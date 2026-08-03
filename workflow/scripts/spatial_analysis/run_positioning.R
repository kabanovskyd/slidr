# set Lib Path
lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path) || lib_path == "") {
  # fallback to default
  .libPaths()[1]
} else {
  .libPaths(lib_path)
}

# load packages
suppressMessages(library(glue)) ; g=glue ; len=length
suppressMessages(library(gridExtra))
suppressMessages(library(magrittr))
suppressMessages(library(Matrix))
suppressMessages(library(jsonlite))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(Seurat))
suppressMessages(library(rlist))
suppressMessages(library(rhdf5))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(qpdf))
suppressMessages(library(qs))
options(warn = -1)
options(future.globals.maxSize = Inf)

# Get the path to the R function
r_func_path <- Sys.getenv('R_FUNC', unset = NA)
if (is.na(r_func_path) || r_func_path == "") {
  # inlined rather than using trouble(): this check runs before any helper file is sourced
  stop(paste0("R_FUNC environment variable is not set, so the helper functions cannot be located.",
              "\nTroubleshooting:",
              "\n \u2022 R_FUNC is exported by the pipeline's spatial-analysis stage; run this script through",
              " `./slidr --bcl <BCL_ID> --spatial-analysis` rather than by hand",
              "\n \u2022 To run it by hand, set it to the scripts' functions/ directory, e.g.",
              " R_FUNC=workflow/scripts/spatial_analysis/functions"), call. = FALSE)
}
source(paste0(r_func_path, '/run_func.R'))
load_matrix_src <- paste0(r_func_path, '/load_matrix.R')
positioning_src <- paste0(r_func_path, '/positioning.R')


# per-sample processing script invoked once per sample by run_spatial.R
# args: RNApath, molecule_info_path, summary_path, SBpath, out_path, ncores, percent_umi_filter
args <- commandArgs(trailingOnly = TRUE)
RNApath <- args[[1]]
molecule_info_path <- args[[2]]
summary_path <- args[[3]]
SBpath <- args[[4]]
out_path <- args[[5]]
ncores <- as.numeric(args[[6]])
# top N percent of beads (by total UMI count) to filter out; defaults to 1 (top 1%) if not provided
percent_umi_filter <- ifelse(length(args) >= 7, as.numeric(args[[7]]), 1)

### Load the RNA ###############################################################
# load_seurat runs normalization, PCA, clustering, and UMAP; attaches cb (barcode without lane suffix)
obj <- load_seurat(RNApath)
# load_intronic reads molecule_info.h5 to add pct.intronic to cell metadata
if (file.exists(molecule_info_path)) { obj %<>% load_intronic(molecule_info_path) }
Misc(obj, "RNA_metadata") = plot_metrics_summary(summary_path, out_path)

### Plots ###############################################################
plot <- UvsI(obj, molecule_info_path)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"RNA.pdf"), 7, 8)))

plot <- plot_umaps(obj)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"UMAP.pdf"), 7, 8)))

### Load the Spatial ###########################################################
# write the cell barcode whitelist (no lane suffix) for load_matrix.R to match against
# then run load_matrix.R as a subprocess to produce matrix.csv.gz and spatial_metadata.json
log_ts("building the spatial matrix")
t_matrix <- Sys.time()
cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))
# shQuote every arg: system2() pastes args unquoted into /bin/sh -c, so metadata-derived paths
# containing shell metacharacters must be escaped to prevent command injection
result <- system2("Rscript", args = shQuote(c(load_matrix_src, SBpath, file.path(out_path, 'cb_whitelist.txt'), out_path, as.character(percent_umi_filter))))
if (result != 0) {
  stop(trouble(
    paste("load_matrix.R failed with code", result),
    "load_matrix.R's own error is above -- it names the real cause",
    paste0("It reads the spatial barcode counts from ", SBpath, "; check that file exists and is complete"),
    "A mismatch between the cell barcodes and the puck barcodes is the usual cause; check the `Puck ID` metadata column names the right puck",
    "Out-of-memory kills also land here on deeply sequenced samples; set `workflow.spatial_downsampling` (e.g. 0.5) to thin the reads",
    paste0("Full log: ", file.path(dirname(dirname(dirname(out_path))), "log", "spatial_analysis.log"))
  ), call. = FALSE)
}
if (!all(file.exists(file.path(out_path, "matrix.csv.gz"), file.path(out_path, "spatial_metadata.json")))) {
  stop(trouble(
    paste0("load_matrix.R reported success but did not write matrix.csv.gz and spatial_metadata.json into ", out_path),
    "Check the output directory is writable and the filesystem has free space",
    "Re-run the spatial analysis from scratch with --spatial-analysis --force"
  ), call. = FALSE)
}
log_ts(sprintf("✓ matrix built — %s", format_duration(difftime(Sys.time(), t_matrix, units = "secs"))))
Misc(obj, "spatial_metadata") <- fromJSON(file.path(out_path, "spatial_metadata.json"))

# run positioning.R as a subprocess; it reads matrix.csv.gz and writes coords.csv
log_ts("positioning cells")
t_position <- Sys.time()
result <- system2("Rscript", args = shQuote(c(positioning_src, file.path(out_path, 'matrix.csv.gz'), out_path, as.character(ncores))))
if (result != 0) {
  stop(trouble(
    paste("positioning.R failed with code", result),
    "positioning.R's own error is above -- it names the real cause",
    "A sample where too few cells carry spatial barcodes cannot be positioned; check the SB UMI counts in the log",
    "Check the puck coordinates are sane -- a puck CSV from a different puck places no cells",
    paste0("Out-of-memory kills also land here; the run used ", ncores, " core(s), and lowering it reduces peak memory"),
    paste0("Full log: ", file.path(dirname(dirname(dirname(out_path))), "log", "spatial_analysis.log"))
  ), call. = FALSE)
}
if (!file.exists(file.path(out_path, "coords.csv"))) {
  stop(trouble(
    paste0("positioning.R reported success but did not write coords.csv into ", out_path),
    "Check the output directory is writable and the filesystem has free space",
    "Re-run the spatial analysis from scratch with --spatial-analysis --force"
  ), call. = FALSE)
}
log_ts(sprintf("✓ cells positioned — %s", format_duration(difftime(Sys.time(), t_position, units = "secs"))))
coords <- read.table(file.path(out_path,"coords.csv"), header=T, sep=",")
# right_join to a full 1:N cb_index frame ensures every called cell is represented,
# even those with no spatial data (they will have NA coordinates)
coords %<>% right_join(data.frame(cb_index=1:len(cb_whitelist)), by = "cb_index") %>% arrange(cb_index)
Misc(obj, "coords") <- coords

# transfer DBSCAN and KDE coordinates onto the Seurat metadata
if (!isTRUE(nrow(coords) == ncol(obj)) || !isTRUE(all(coords$cb_index == obj$cb_index))) {
  stop(trouble(
    paste0("the positioned coordinates do not line up with the Seurat object (coords rows: ", nrow(coords),
           ", object cells: ", ncol(obj), ")"),
    "The coordinate table is joined onto the cells by cb_index, so a mismatch would silently assign cells the wrong positions",
    paste0("This usually means ", out_path, " holds coords.csv from an earlier run against a different count matrix"),
    "Re-run the spatial analysis from scratch with --spatial-analysis --force",
    "If it persists after a clean re-run, please report it at https://github.com/kabanovskyd/slidr/issues"
  ), call. = FALSE)
}
obj$x_um <- coords$x_um
obj$y_um <- coords$y_um
obj$x_um_dbscan <- coords$x_um_dbscan
obj$y_um_dbscan <- coords$y_um_dbscan

# "dbscan" reduction: raw DBSCAN centroid for all placed cells (cluster==1 only)
emb = coords %>% select(x_um_dbscan, y_um_dbscan)
log_detail(sprintf("%d of %d cells placed by DBSCAN",
                   sum(!is.na(emb$x_um_dbscan)), nrow(emb)))
colnames(emb) = c("d_1","d_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["dbscan"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "d_", assay = "RNA")

# "kde" reduction: KDE peak location; masked to NA when ratio > 1/3 (ambiguous placement)
emb = coords %>% mutate(across(everything(), ~ifelse(ratio > 1/3, NA, .))) %>% select(x_um_kde, y_um_kde)
colnames(emb) = c("k_1","k_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["kde"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "k_", assay = "RNA")

# "spatial" reduction: the final placement (KDE-filtered DBSCAN — NA when ratio >= 1/3)
# this is the recommended reduction for downstream spatial analysis
emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

# save Seurat object in qs format (faster and smaller than RDS for large objects)
qsave(obj, file.path(out_path,"seurat.qs"))

# Plots DBSCAN
plot <- plot_clusters(obj, reduction="dbscan")
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"DimPlotDBSCAN.pdf"), 7, 8)))

# Plots KDE
plot <- plot_clusters(obj, reduction="kde")
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"DimPlotKDE.pdf"), 7, 8)))

# Plots DimPlot
plot <- plot_clusters(obj, reduction="spatial")
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path,"DimPlot.pdf"), 7, 8)))

# Plots RNA vs SB
plot <- plot_RNAvsSB(obj)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "RNAvsSB.pdf"), 7, 8)))

# merge all per-stage PDFs into a single summary.pdf in a logical reading order:
# (RNA QC → spatial QC → UMAP → spatial maps → diagnostic plots)
# only include files that were actually generated (some are optional depending on data availability)
plotlist <- c(c("SB.pdf","beadplot.pdf","SBmetrics.pdf"),
              c("DBSCAN.pdf","KDE.pdf","DBSCANvsKDE.pdf","beadplots.pdf"),
              c("RNAmetrics.pdf","RNA.pdf","UMAP.pdf","DimPlot.pdf","DimPlotDBSCAN.pdf","DimPlotKDE.pdf","RNAvsSB.pdf"),
              c("beadfilter.pdf"))
plotorder <- c(8, 9, 10, 1, 2, 15, 4, 5, 6, 11, 12, 13, 14, 3, 7)
suppressMessages({
  pdfs <- file.path(out_path, plotlist[plotorder])
  pdfs %<>% keep(file.exists)
  qpdf::pdf_combine(input=pdfs, output=file.path(out_path,"summary.pdf"))
  file.remove(pdfs)
})
