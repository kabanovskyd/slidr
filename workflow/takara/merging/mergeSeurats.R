rm(list=ls())
library(Seurat)
library(Matrix)
library(sceasy)

options(future.globals.maxSize = 5000 * 1024^2)

# locate this script's own directory so the shared helpers file can be found regardless of cwd
.script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.script_file_arg) != 1) {
  stop(paste0("unable to determine script directory: mergeSeurats.R must be run via 'Rscript --vanilla <path>'",
       "\nTroubleshooting:",
       "\n \u2022 The script locates its shared helpers relative to its own path, which it reads from the --file= argument",
       "\n \u2022 Sourcing it from an interactive session or piping it on stdin removes that argument",
       "\n \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`"), call.=FALSE)
}
source(file.path(dirname(normalizePath(sub("^--file=", "", .script_file_arg))), "..", "shared_seurat_utils.R"))

############ FUNCTIONS #############
rotate_spatial<-function(coords){
  data.frame(x=-coords[,2], y=coords[,1], row.names=rownames(coords), stringsAsFactors=FALSE)
}

# run analysis
analyze <- function(obj, res=0.2, n.epochs=NULL) {
  obj <- Seurat::SCTransform(obj, assay="RNA", ncells=1000, verbose=T, conserve.memory=T)
  obj <- Seurat::RunPCA(obj, verbose=T)
  npc <- ncol(Loadings(obj, reduction="pca"))
  ndim <- min(30, npc)
  # seq_len(ndim) instead of 1:ndim so a degenerate PCA (npc==0) yields an empty vector rather
  # than the c(1,0) reverse-range footgun
  obj <- Seurat::RunUMAP(obj, dims=seq_len(ndim))
  obj <- Seurat::FindNeighbors(obj, dims=seq_len(ndim))
  obj <- Seurat::FindClusters(obj, resolution=res)
  obj  
}

# add metadata
add.meta <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  obj[["logumi"]] <- log10(obj$nCount_RNA+1)
  return(obj)
}

# organize per cluster markers
gen_markers_bycluster<-function(variable_features_cluster){
  markers_bycluster<-lapply(unique(variable_features_cluster$cluster), function(x){
    genes<-variable_features_cluster[variable_features_cluster$cluster==x,]
    genes.top<-genes[genes$p_val_adj<0.05,]
    genes.top<-genes.top[order(genes.top$avg_log2FC, decreasing = T),]
  });names(markers_bycluster)<-unique(variable_features_cluster$cluster)
  return(markers_bycluster)
}

############ INPUT CHECKS ############
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop(paste0("insufficient command-line arguments provided.",
       "\n Usage: Rscript --vanilla mergeSeurats.R <input_folder_path> <pattern> <output_prefix> <outdir> <profile_conda_path>",
       "\nTroubleshooting:",
       "\n \u2022 This script is normally invoked by the Flex spatial-analysis stage, not by hand",
       "\n \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
       "\n \u2022 Received ", length(args), " argument(s): ", paste(args, collapse=" ")), call.=FALSE)
}else{
  input_folder_path <- args[1]
  pattern <- args[2]
  output_prefix <- args[3]
  outdir <- args[4]
  profile_conda_path = args[5]
}

intermediates_dir <- file.path(outdir, "intermediates")
# use full.names instead of setwd(input_folder_path) -- on.exit() at top-level script scope
# (no enclosing function) never fires when run via Rscript, so a setwd() here would leave the
# process's cwd permanently changed; using full paths avoids touching global state at all
seurat_files <- list.files(path = input_folder_path, pattern = pattern, full.names = TRUE)
if (length(seurat_files) == 0) {
  stop(paste0("no files matching pattern '", pattern, "' found in ", input_folder_path,
       "\nTroubleshooting:",
       "\n \u2022 This step merges the per-partition Seurat objects produced by the Flex spatial-profiling step",
       "\n \u2022 Nothing to merge means none of this sample's partitions completed profiling",
       "\n \u2022 Check the per-partition logs in the run's takara_pipeline.log for the underlying error",
       "\n \u2022 List what is actually there: `ls -l ", input_folder_path, "`"), call.=FALSE)
}

library(reticulate)
reticulate::use_condaenv(profile_conda_path, required = TRUE)

message("########## Load each Seurat object into a list ##########")
seurat_list <- list()
for (file in seurat_files) {
  seurat_obj <- readRDS(file)
  file_prefix <- gsub("_+$", "", gsub(pattern,"",basename(file)))
  message(file_prefix, " has ", length(Cells(seurat_obj)), " nuclei")
  seurat_list[[file_prefix]] <- seurat_obj
}

message("################## Merge Seurat objects #################")
# guard the single-object case: seurat_list[-1] would be an empty list, and merge(x, y=list())
# is undefined; use the lone object directly (item 69 already handles the zero-file case upstream)
if (length(seurat_list) < 2) {
  message("Only one Seurat object found; skipping merge and using it directly")
  merged_sobj <- seurat_list[[1]]
} else {
  merged_sobj <- merge(seurat_list[[1]], y = seurat_list[-1], merge.dr = TRUE, add.cell.ids = names(seurat_list))
}
message("Merged seurat object has ", length(Cells(merged_sobj)), " nuclei")
# Remove spatial info stored in "@images$slices*" except for slice1
merged_sobj@images[-1]=NULL

message("#### Load merged spatial coord into '@images$slice1' ####")
coords <- rotate_spatial(as.data.frame(merged_sobj[["SPATIAL"]]@cell.embeddings))
merged_sobj@images$slice1 = new(Class = "SlideSeq", assay = "Spatial", key = "slice1_", coordinates = coords)
s_dims <- dim(merged_sobj@images$slice1@coordinates)
message("dimension of @images$slice1: ", s_dims[1], " ", s_dims[2])

message("################### Rerun analysis ######################")
merged_sobj <- analyze(merged_sobj, res=0.2)
variable_features_cluster <- FindAllMarkers(merged_sobj, assay = "SCT", only.pos =T)
markers_bycluster <- gen_markers_bycluster(variable_features_cluster)

message("# Generate Positioned and Conf Positioned Seurat objects #")
merged_sobj_pos <- subset(merged_sobj, cells=Cells(merged_sobj)[merged_sobj$number_clusters>0])
message("Number of postioned nuclei: ", length(Cells(merged_sobj_pos)))
if (length(Cells(merged_sobj_pos)) == 0) {
  stop(paste0("no nuclei passed QC filtering (number_clusters>0) for merged sample ", output_prefix,
       "\nTroubleshooting:",
       "\n \u2022 No nucleus in any partition of this sample was assigned a spatial location",
       "\n \u2022 Check the `Puck ID` metadata column names the puck this sample was actually run on -- the wrong puck places no nuclei",
       "\n \u2022 Check the spatial-barcode library was sequenced deeply enough; the read and UMI counts are in the run's takara_pipeline.log",
       "\n \u2022 Check the `Flex Probe Barcode IDs` metadata column matches the barcodes actually used, so reads are not all discarded"), call.=FALSE)
}

merged_sobj_pos_conf <- subset(merged_sobj_pos, cells=Cells(merged_sobj_pos)[merged_sobj_pos$number_clusters==1])
message("Number of confidently positioned nuclei: ", (length(Cells(merged_sobj_pos_conf))))
if (length(Cells(merged_sobj_pos_conf)) == 0) {
  stop(paste0("no nuclei passed QC filtering (number_clusters==1) for merged sample ", output_prefix,
       "\nTroubleshooting:",
       "\n \u2022 Nuclei were placed, but none of them at a single unambiguous location",
       "\n \u2022 A `number_clusters==1` filter keeps only unambiguously placed nuclei, so a sample with diffuse spatial signal can pass the previous check and fail this one",
       "\n \u2022 Check the `Puck ID` metadata column names the puck this sample was actually run on -- the wrong puck places no nuclei",
       "\n \u2022 Check the spatial-barcode library was sequenced deeply enough; the read and UMI counts are in the run's takara_pipeline.log",
       "\n \u2022 Check the `Flex Probe Barcode IDs` metadata column matches the barcodes actually used, so reads are not all discarded"), call.=FALSE)
}

## re-find spatially variable features on conf positioned nuclei
merged_sobj_pos_conf <- FindSpatiallyVariableFeatures(merged_sobj_pos_conf, assay="SCT", slot = "scale.data",
						      features = head(VariableFeatures(merged_sobj_pos_conf), 200),
						      selection.method = "moransi", x.cuts = 100, y.cuts = 100,
						      verbose = TRUE, nfeatures=200)
variable_features_moransi <- organize_spatial_features_moransi(merged_sobj_pos_conf)

message("############### Save merged Seurat objects ###############")
saveRDS(merged_sobj, file=file.path(intermediates_dir, paste0(output_prefix, "_seurat_spatial_merged.rds")))
saveRDS(merged_sobj_pos, file=file.path(intermediates_dir, paste0(output_prefix, "_Positioned_seurat_spatial_merged.rds")))
saveRDS(merged_sobj_pos_conf, file=file.path(outdir, paste0(output_prefix, "_ConfPositioned_seurat_spatial_merged.rds")))

message("################ Save merged H5AD objects ################")
### adjust spatial coord before converting to h5ad object
merged_sobj <- flip_spatial(merged_sobj)
merged_sobj_pos <- flip_spatial(merged_sobj_pos)
merged_sobj_pos_conf <- flip_spatial(merged_sobj_pos_conf)
sceasy::convertFormat(merged_sobj, from="seurat", to="anndata", outFile=file.path(intermediates_dir, paste0(output_prefix, "_anndata_merged.h5ad")))
sceasy::convertFormat(merged_sobj_pos, from="seurat", to="anndata", outFile=file.path(intermediates_dir, paste0(output_prefix, "_Positioned_anndata_merged.h5ad")))
sceasy::convertFormat(merged_sobj_pos_conf, from="seurat", to="anndata", outFile=file.path(outdir, paste0(output_prefix, "_ConfPositioned_anndata_merged.h5ad")))

message("################## Save sparse matrices #################")
merged_coords <- as.data.frame(merged_sobj[["SPATIAL"]]@cell.embeddings)
output_sparse_matrix(output_prefix, merged_sobj, merged_coords, "_merged", intermediates_dir)
output_sparse_matrix(output_prefix, merged_sobj_pos, merged_coords, "_PositionedNuclei_merged", intermediates_dir)
output_sparse_matrix(output_prefix, merged_sobj_pos_conf, merged_coords, "_ConfPositionedNuclei_merged", outdir)

message("################## Save variable genes ##################")
write.table(variable_features_moransi, file.path(outdir, paste0(output_prefix, "_variable_features_spatial_moransi_merged.txt")), quote = F, sep="\t")
write.csv(variable_features_cluster, file.path(outdir, paste0(output_prefix, "_variable_features_clusters_merged.csv")), quote=F)
