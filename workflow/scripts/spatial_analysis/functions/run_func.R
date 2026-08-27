#########################################################################################
################################ RUN JOBS FUNCTIONS ##################################
#########################################################################################

source(paste0(r_func_path, '/plot_utils_func.R'))


# gex_summary_format converts the 10X multi summary.csv column names (GEX. prefixed)
# to the standard Cellranger metrics_summary.csv column names expected by downstream plotting
gex_summary_format <- function(summary_path) {
  df <- read.csv(summary_path)
  column_mapping <- c(
    'Estimated.number.of.cells' = 'Estimated Number of Cells',
    'GEX.Mean.raw.reads.per.cell' = 'Mean Reads per Cell',
    'GEX.Median.genes.per.cell' = 'Median Genes per Cell',
    'GEX.Sequenced.read.pairs' = 'Number of Reads',
    'GEX.Valid.barcodes' = 'Valid Barcodes',
    'GEX.Percent.duplicates' = 'Sequencing Saturation',
    'GEX.Q30.bases.in.barcode' = 'Q30 Bases in Barcode',
    'GEX.Q30.bases.in.read.2' = 'Q30 Bases in RNA Read',
    'GEX.Q30.bases.in.UMI' = 'Q30 Bases in UMI',
    'GEX.Reads.mapped.to.genome' = 'Reads Mapped to Genome',
    'GEX.Reads.mapped.confidently.to.genome' = 'Reads Mapped Confidently to Genome',
    'GEX.Reads.mapped.confidently.to.intergenic.regions' = 'Reads Mapped Confidently to Intergenic Regions',
    'GEX.Reads.mapped.confidently.to.intronic.regions' = 'Reads Mapped Confidently to Intronic Regions',
    'GEX.Reads.mapped.confidently.to.exonic.regions' = 'Reads Mapped Confidently to Exonic Regions',
    'GEX.Reads.mapped.confidently.to.transcriptome' = 'Reads Mapped Confidently to Transcriptome',
    'GEX.Reads.mapped.antisense.to.gene' = 'Reads Mapped Antisense to Gene',
    'GEX.Fraction.of.transcriptomic.reads.in.cells' = 'Fraction Reads in Cells',
    'GEX.Total.genes.detected' = 'Total Genes Detected',
    'GEX.Median.UMI.counts.per.cell' = 'Median UMI Counts per Cell'
  )
  extract_info_columns <- names(column_mapping)
  df_extract <- df %>% select(all_of(extract_info_columns))
  colnames(df_extract) <- column_mapping
  
  percentage_columns <- c(
    'Valid Barcodes',
    'Sequencing Saturation',
    'Q30 Bases in Barcode',
    'Q30 Bases in RNA Read',
    'Q30 Bases in UMI',
    'Reads Mapped to Genome',
    'Reads Mapped Confidently to Genome',
    'Reads Mapped Confidently to Intergenic Regions',
    'Reads Mapped Confidently to Intronic Regions',
    'Reads Mapped Confidently to Exonic Regions',
    'Reads Mapped Confidently to Transcriptome',
    'Reads Mapped Antisense to Gene',
    'Fraction Reads in Cells'
  )
  df_extract <- df_extract %>%
    mutate(across(all_of(percentage_columns), ~ paste0(. * 100, "%")))
  
  # sequencing saturation = 1 - (unique reads / total reads); the original column is percent duplicates
  # so unique reads = total * (1 - pct_dup) and saturation = 1 - unique/total = pct_dup
  tr <- df$GEX.Sequenced.read.pairs
  ur <- tr * (1 - df$GEX.Percent.duplicates)
  df_extract$`Sequencing Saturation` <- paste0((1 - (ur / tr)) * 100, "%")
  
  # save name
  save_path <- gsub('summary.csv', 'metrics_summary.csv', summary_path)
  write.csv(df_extract, save_path, row.names = FALSE)
}


check_args_and_get_paths <- function(args) {
  if (length(args) < 4) {
    stop(trouble(
      "Usage: Rscript run_spatial.R BCLname [Sample1, Sample2] <cellbender> <downsampling_rate> [ncores] [percent_umi_filter]",
      "This script is normally invoked by the pipeline's spatial-analysis stage, not by hand",
      "Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
      paste("Got", length(args), "argument(s):", paste(args, collapse = " "))
    ), call. = FALSE)
  }
  # DATA_PATH is set by the pipeline to the base output directory containing all BCL run folders
  base_path <- Sys.getenv('DATA_PATH', unset = NA)
  if (is.na(base_path) || base_path == "") {
    stop(trouble(
      "DATA_PATH environment variable is not set, so the run's output directory cannot be located.",
      "DATA_PATH is exported by the pipeline's spatial-analysis stage; run this script through `./slidr --bcl <BCL_ID> --spatial-analysis` rather than by hand",
      "To run it by hand, set it to the directory that contains the per-BCL output folders (the `output_path` from your config file)"
    ), call. = FALSE)
  }

  # SAMPLEnames is passed as a bracket-delimited string "[Sample1, Sample2]" to survive shell quoting
  BCLname <- args[1]
  SAMPLEnames <- args[2]
  cellbender <- as.logical(args[3])
  downsampling_rate <- args[4]
  ncores <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
  # top N percent of beads (by total UMI count) to filter out; defaults to 1 (top 1%) if not provided
  percent_umi_filter <- ifelse(length(args) >= 6, as.numeric(args[6]), 1)

  if (!is.null(SAMPLEnames)) {
    if (!grepl("^\\[.*\\]$", SAMPLEnames)) {
      stop(trouble(
        paste0("Invalid format for sample names: '", SAMPLEnames, "'"),
        "The sample list must be bracket-delimited, e.g. '[Sample1, Sample2]' -- the brackets are what survives shell quoting",
        "The pipeline builds this argument itself, so seeing it here usually means the script was invoked by hand",
        "Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`"
      ), call. = FALSE)
    }
    SAMPLEnames <- gsub("\\[|\\]", "", SAMPLEnames)
    SAMPLEnames <- strsplit(SAMPLEnames, ",\\s*")[[1]]
  }

  RNApath = paste0(base_path, "/", BCLname, "/output/count")
  SBpath = paste0(base_path, "/", BCLname, "/output/spatial_barcodes")
  if (cellbender) {
    cellbender_path <- paste0(base_path, "/", BCLname, "/output/cellbender")
  }

  if (!file.exists(RNApath)) {
    stop(trouble(
      paste0("no cellranger count output directory found for BCL '", BCLname, "' at ", RNApath),
      "Run the count stage first with --count (or the whole pipeline with --run-all)",
      "Check --bcl names the run whose outputs you want to analyse",
      paste0("Check `paths.output_path` in your config file resolves to ", base_path)
    ), call. = FALSE)
  }
  subfolders <- list.files(RNApath, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  names(subfolders) <- sapply(subfolders, basename)
  # if no sample names provided, process all subfolders found in the count output directory
  if (is.null(SAMPLEnames) | length(SAMPLEnames) == 0) {SAMPLEnames <- basename(subfolders)}

  files_df <- data.frame(name = character(),
                         rna = character(),
                         molecule = character(),
                         summary = character(),
                         spatial = character(),
                         output = character(),
                         ncores = numeric(),
                         percent_umi_filter = numeric(),
                         stringsAsFactors = FALSE)

  sample_rows <- vector("list", length(SAMPLEnames))

  for (i in seq_along(SAMPLEnames)) {
    name <- SAMPLEnames[i]
    if (cellbender) {
      rna_file <- paste0(cellbender_path, "/", name, "/cellbender_output_filtered.h5")
    } else {
      rna_file <- paste0(RNApath, "/", name, "/filtered_feature_bc_matrix.h5")
    }
    mol_file <- paste0(RNApath, "/", name, "/molecule_info.h5")
    summary_file <- paste0(RNApath, "/", name, "/metrics_summary.csv")
    # 10X multi outputs use summary.csv and gex_molecule_info.h5 instead of the standard names
    gex_summary_file <- paste0(RNApath, "/", name, "/summary.csv")
    gex_mol_file <- paste0(RNApath, "/", name, "/gex_molecule_info.h5")
    if (downsampling_rate == "") {
      sb_file <- paste0(SBpath, "/", name, "/SBcounts.h5")
    } else {
      sb_file <- paste0(SBpath, "/", name, "/SBcounts_downsampled_", downsampling_rate, ".h5")
    }

    rna_path <- ifelse(file.exists(rna_file), rna_file, NA)
    mol_path <- ifelse(file.exists(mol_file), mol_file, 
                       ifelse(file.exists(gex_mol_file), gex_mol_file, NA))
    summary_path <- ifelse(file.exists(summary_file), summary_file, 
                          ifelse(file.exists(gex_summary_file), gex_summary_file, NA))
    sb_path <- ifelse(file.exists(sb_file), sb_file, NA)

    missing_files <- c()
    if (is.na(rna_path)) {
      if (cellbender) {
        missing_files <- c(missing_files, "/cellbender_output_filtered.h5")
      } else {
        missing_files <- c(missing_files, "filtered_feature_bc_matrix.h5")
      }
    }
    if (is.na(mol_path)) {
      missing_files <- c(missing_files, "molecule_info.h5 or gex_molecule_info.h5")
    }
    if (is.na(summary_path)) {
      missing_files <- c(missing_files, "metrics_summary.csv")
    }
    if (is.na(sb_path)) {
      # the name that was looked for, not the literal "SBcounts.h5": under a downsampling rate the
      # file being asked for is SBcounts_downsampled_<rate>.h5, and reporting the undownsampled name
      # sends the reader to a file that is sitting right there
      missing_files <- c(missing_files, basename(sb_file))
    }
    if (length(missing_files) > 0) {
      # name the stage that produces each missing file, so the message says which step to re-run
      hints <- c()
      if (any(grepl("cellbender", missing_files))) {
        hints <- c(hints, "Run the cellbender stage first with --cellbender, or pass --no-cellbender to analyse the raw cellranger counts instead")
      }
      if (any(grepl("filtered_feature_bc_matrix|molecule_info|metrics_summary", missing_files))) {
        hints <- c(hints, paste0("Run the count stage first with --count; its outputs belong in ", RNApath, "/", name, "/"))
      }
      if (any(grepl("SBcounts", missing_files))) {
        hints <- c(hints, paste0("Run the spatial-barcode counting stage first with --spatial-count; its output belongs in ", SBpath, "/", name, "/"))
        if (downsampling_rate != "") {
          hints <- c(hints, paste0("This run expects the downsampled file SBcounts_downsampled_", downsampling_rate,
                                   ".h5 because `workflow.spatial_downsampling` is set; clear that field to use the full SBcounts.h5"))
        }
      }
      stop(trouble(
        paste0("required input files are missing for sample '", name, "': ", paste(missing_files, collapse = ", ")),
        c(hints,
          "Check the `Sample Name` metadata column matches the output directory names exactly",
          "An earlier stage that crashed part-way leaves an incomplete directory; re-run it with --force")
      ), call. = FALSE)
    }

    log_ts(paste("sample", name))
    log_detail(c(paste("rna      →", rna_path),
                 paste("molecule →", mol_path),
                 paste("spatial  →", sb_path)))

    out_path <- gsub("spatial_barcodes", "spatial_analysis", gsub(basename(sb_path), '', sb_path))

    if (!file.exists(out_path)) {dir.create(out_path, recursive = TRUE)}
    if (summary_path == gex_summary_file) {
      gex_summary_format(summary_path)
    }
    sample_rows[[i]] <- data.frame(name = name,
               rna = rna_path,
               molecule = mol_path,
               summary = summary_path,
               spatial = sb_path,
               output = out_path,
               ncores = ncores,
               percent_umi_filter = percent_umi_filter)
  }
  files_df <- do.call(rbind, c(list(files_df), sample_rows))
  return(files_df)
}




# ReadCB_h5 reads a CellBender output HDF5 file, which uses a different internal structure
# than the 10X filtered_feature_bc_matrix.h5 (CellBender stores data under a "matrix" group)
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop(trouble(
      "the hdf5r package is required to read CellBender .h5 files, but is not installed",
      "The pipeline installs it into the `r_env` conda environment; recreate that environment with `mamba env remove -n r_env` then re-run",
      "Or install it into the active library by hand: install.packages('hdf5r')",
      paste0("Active R library path(s): ", paste(.libPaths(), collapse = ", "))
    ), call. = FALSE)
  }
  if (!file.exists(filename)) {
    stop(trouble(
      paste0("CellBender output file not found: ", filename),
      "Run the cellbender stage first with --cellbender (or the whole pipeline with --run-all)",
      "Or pass --no-cellbender to analyse the raw cellranger count matrices instead",
      "A cellbender run that crashed part-way leaves no filtered .h5; re-run it with --cellbender --force"
    ), call. = FALSE)
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  
  if ("matrix" %in% genomes) {
    genome <- "matrix"
  } else {
    stop(trouble(
      paste0("no 'matrix' group in the CellBender file ", filename,
             " (found: ", paste(genomes, collapse = ", "), ")"),
      "CellBender writes its counts under a top-level 'matrix' group; a file without one is not a CellBender output",
      "Check this is cellbender_output_filtered.h5 and not a cellranger matrix or a truncated file",
      "Re-run the cellbender stage with --cellbender --force to regenerate it",
      "Or pass --no-cellbender to analyse the raw cellranger count matrices instead"
    ), call. = FALSE)
  }
  if (hdf5r::existsGroup(infile, "matrix")) {
    if (use.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }
  else {
    if (use.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  counts <- infile[[paste0(genome, "/data")]]
  indices <- infile[[paste0(genome, "/indices")]]
  indptr <- infile[[paste0(genome, "/indptr")]]
  shp <- infile[[paste0(genome, "/shape")]]
  features <- infile[[paste0(genome, "/", feature_slot)]][]
  barcodes <- infile[[paste0(genome, "/barcodes")]]
  sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                             x = as.numeric(x = counts[]), dims = shp[], repr = "T")
  if (unique.features) {features <- make.unique(names = features)}
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as.sparse(x = sparse.mat)
  if (infile$exists(name = paste0(genome, "/features"))) {
    types <- infile[[paste0(genome, "/features/feature_type")]][]
    types.unique <- unique(x = types)
    if (length(x = types.unique) > 1) {
      log_detail(paste0("genome ", genome, " has multiple modalities — returning a list of matrices"))
      sparse.mat <- sapply(X = types.unique, FUN = function(x) {
        return(sparse.mat[which(x = types == x), ])
      }, simplify = FALSE, USE.NAMES = TRUE)
    }
  }
  infile$close_all()
  return(sparse.mat)
}



load_seurat <- function(RNApath) {
  log_ts("running Seurat processing (normalize → PCA → cluster → UMAP)")
  if(str_ends(RNApath, "filtered_feature_bc_matrix.h5")) {
    obj_data <- Seurat::Read10X_h5(RNApath)
    # multiome libraries contain both Gene Expression and Peaks assays
    if (all(c("Gene Expression", "Peaks") %in% names(obj_data))) {
      obj <- CreateSeuratObject(counts = obj_data$`Gene Expression`)
      obj[["Peaks"]] <- CreateAssayObject(counts = obj_data$Peaks)
    } else {
      obj <- CreateSeuratObject(obj_data)
    }
  } else {
    obj <- ReadCB_h5(RNApath) %>% CreateSeuratObject()
  }
  # strip the lane suffix (e.g. "-1") from Cell Ranger barcodes to get the raw barcode sequence
  # cb_index is a 1-based integer that links Seurat cells to the SB count matrix cb_index
  obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .)) %T>% {stopifnot(!duplicated(.))}
  obj[["cb_index"]] <- 1:ncol(obj)
  obj[["logumi"]] <- log10(obj$nCount_RNA+1)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  
  # Add tech metadata
  if (all(c("Peaks", "RNA") %in% names(obj@assays))) {
    Misc(obj, "tech") <- "Multiome"
  } else if (identical(names(obj@assays), "RNA")) {
    Misc(obj, "tech") <- "RNA"
  } else {
    stop(trouble(
      paste0("the Seurat object has neither the expected 'RNA' assay alone nor an 'RNA'+'Peaks' pair (found: ",
             paste(names(obj@assays), collapse = ", "), ")"),
      paste0("Check the count matrix at ", RNApath, " is a gene-expression (or multiome) matrix"),
      "A feature-barcoding or antibody-capture library produces different assay names and is not supported by this stage",
      "Re-run the count stage with --count --force if the matrix looks wrong"
    ), call. = FALSE)
  }
  
  # PCA, Cluster, and UMAP
  obj %<>% Seurat::NormalizeData(verbose = FALSE) %>%
    Seurat::FindVariableFeatures(verbose = FALSE) %>%
    Seurat::ScaleData(verbose = FALSE) %>%
    Seurat::RunPCA(npcs=50, verbose = FALSE) 
  gc()
  
  obj %<>% Seurat::FindNeighbors(dims=1:30, verbose = FALSE) %>%
    Seurat::FindClusters(resolution=round(ncol(obj) / 20000, 2), verbose = FALSE)
  suppressMessages(suppressWarnings(
    obj <- Seurat::RunUMAP(obj, dims=1:30, verbose = FALSE, n.epochs=NULL)
  ))
  
  log_detail(paste0(ncol(obj), " cells"))
  Misc(obj, "RNApath") <- RNApath

  return(obj)
}



load_intronic <- function(obj, molecule_info_path) {
  library(rhdf5)
  fetch <- function(x){return(h5read(molecule_info_path, x))}
  # umi_type==0 indicates an intronic UMI in the Cell Ranger molecule_info.h5 schema
  # barcode_idx is 0-based; +1 converts to 1-based indexing into the barcodes array
  if ("barcode_idx" %in% h5ls(molecule_info_path)$name) {
    barcodes = fetch("barcodes")
    info = data.frame(barcode=fetch("barcode_idx")+1, umi_type=fetch("umi_type"))
    info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
    obj$pct.intronic = info$pct.intronic[match(obj$cb, barcodes[info$barcode])] * 100
  } else {
    log_ts("[WARNING]: no intronic information found in the molecule_info file")
  }
  return(obj)
}



# Plot metrics summary
plot_metrics_summary <- function(summary_path, out_path) {
  if (nchar(summary_path)==0 || !file.exists(summary_path)) {
    make.pdf(gdraw("No metrics_summary.csv found"), file.path(out_path,"RNAmetrics.pdf"), 7, 8)
    return(c())
  }
  plotdf = read.table(summary_path, header=F, sep=",", comment.char="")
  plotdf %<>% t
  rownames(plotdf) = NULL
  
  plot = plot_grid(ggdraw()+draw_label(""),
                   ggdraw()+draw_label(g("10X Metrics Summary")),
                   plot.tab(plotdf),
                   ggdraw()+draw_label(""),
                   ncol=1, rel_heights=c(0.1,0.1,0.7,0.2))
  make.pdf(plot, file.path(out_path, "RNAmetrics.pdf"), 7, 8)
  
  return(setNames(plotdf[,2], plotdf[,1]))
}



# Plot RNA curves
UvsI <- function(obj, molecule_info_path) {
  if (!file.exists(molecule_info_path) || nchar(molecule_info_path) == 0) {
    plot <- gdraw("No molecule_info.h5 found")
    return(plot)
  }
  if (!"barcode_idx" %in% h5ls(molecule_info_path)$name) {
    plot <- gdraw("Unrecognized molecule_info.h5")
    return(plot)
  }
  
  fetch <- function(x){return(h5read(molecule_info_path,x))}
  barcodes = fetch("barcodes")
  molecule_info = data.frame(barcode=fetch("barcode_idx"),
                             umi_type=fetch("umi_type"),
                             reads=fetch("count"))
  
  # Panel 3: downsampling curve
  tab = table(molecule_info$reads)
  downsampling = map_int(seq(0,1,0.05), function(p){
    map2_int(tab, as.numeric(names(tab)), function(v, k){
      length(unique(floor(sample(0:(k*v-1), round(k*v*p), replace=F)/k)))
    }) %>% sum
  })
  plotdf = data.frame(
    x = seq(0,1,0.05)*sum(molecule_info$reads)/1000/1000,
    y = downsampling/1000/1000
  )
  p0 = ggplot(plotdf, aes(x=x,y=y))+geom_line()+theme_bw()+xlab("Millions of reads")+ylab("Millions of filtered UMIs")+ggtitle("RNA Downsampling curve")
  
  df = molecule_info %>% group_by(barcode) %>% summarize(umi=n(), pct.intronic=sum(umi_type==0)/umi) %>% 
    ungroup %>% arrange(desc(umi)) %>% mutate(logumi=log10(umi))
  
  # Panel 2 and 4: intronic density
  if (!is.null(df$pct.intronic) && !all(df$pct.intronic==0)) {
    ct = 500
    if (any(df$umi>=ct)) {
      p1 = df %>% filter(umi>=ct) %>% ggplot(aes(x = logumi, y = pct.intronic)) + 
        geom_bin2d(bins=100) +
        scale_fill_viridis(trans="log", option="A", name="density") + 
        theme_minimal() +
        labs(title = g("Intronic vs. UMI droplets (>{ct} umi)"), x = "logumi", y = "%intronic") + NoLegend()
      
      max_density_x = density(filter(df,umi>=ct,pct.intronic>1/3)$pct.intronic) %>% {.$x[which.max(.$y)]}
      p2 = df %>% filter(umi>=ct) %>% ggplot(aes(x = pct.intronic)) + geom_density() + 
        theme_minimal() + labs(title = g("Intronic density (>{ct} umi)"), x = "%intronic", y = "Density") + 
        geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
        annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
    } else {
      p1 = gdraw(g("No cells with {ct}+ UMI"))
      p2 = gdraw(g("No cells with {ct}+ UMI"))
    }
  } else {
    p1 = gdraw(g("No intronic information"))
    p2 = gdraw(g("No intronic information"))
  }
  
  # Panel 1: cell barcode knee plot
  df %<>% mutate(index=seq_len(nrow(df)), called=barcodes[barcode+1] %in% obj$cb)
  p3 = ggplot(df,aes(x=index,y=umi,col=called))+geom_line()+theme_bw()+scale_x_log10()+scale_y_log10()+
    ggtitle("Barcode rank plot")+xlab("Cell barcodes")+ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
  
  plot = plot_grid(p3, p1, p0, p2, ncol=2)
  return(plot)
}


# Plot UMAP + metrics
plot_umaps <- function(obj) {
  mytheme <- function(){theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="top", legend.justification="center", legend.key.width=unit(2, "lines"))}
  umap <- DimPlot(obj,label=T) + ggtitle(g("UMAP")) + mytheme() + NoLegend()  + coord_fixed(ratio=1)
  logumi <- VlnPlot(obj,"logumi",alpha=0) + mytheme() + NoLegend() 
  mt <- FeaturePlot(obj,"percent.mt") + ggtitle("%MT") + mytheme() + coord_fixed(ratio=1) + 
    annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$percent.mt), 2)}%\nMean: {round(mean(obj$percent.mt), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  if ("pct.intronic" %in% names(obj@meta.data)) {
    intronic <- FeaturePlot(obj,"pct.intronic") + ggtitle("%Intronic") + mytheme() + coord_fixed(ratio=1) +
      annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$pct.intronic), 2)}%\nMean: {round(mean(obj$pct.intronic), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  } else {
    intronic <- gdraw("No intronic information")
  }
  
  plot <- plot_grid(umap, logumi, mt, intronic, ncol=2)
  return(plot)
}


# Create DimPlot
plot_clusters <- function(obj, reduction) {
  # infer the aspect ratio of the puck layout to choose the number of columns for the per-cluster grid
  npucks = (max(obj$x_um,na.rm=T)-min(obj$x_um,na.rm=T))/(max(obj$y_um,na.rm=T)-min(obj$y_um,na.rm=T))
  nclusters = len(unique(obj$seurat_clusters))
  # when no cells are placed (all x_um/y_um are NA), max()/min() with na.rm return -Inf/Inf and
  # npucks becomes NaN, making ncols NaN -> DimPlot(ncol=NaN) errors; fall back to a single column
  if (!is.finite(npucks) || npucks <= 0) {
    ncols = 1
  } else {
    ncols = round(sqrt(npucks*nclusters/2)/npucks*2)
  }
  if (!is.finite(ncols) || ncols < 1) ncols = 1
  
  m = obj@reductions[[reduction]]@cell.embeddings %>% {!is.na(.[,1]) & !is.na(.[,2])}
  title = g("%placed: {round(sum(m)/len(m)*100,2)} ({sum(m)}/{len(m)}) [{reduction}]")
  p1 = DimPlot(obj, reduction=reduction) + coord_fixed(ratio=1) +
    ggtitle(title) + NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)") + 
    theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=10)) +
    theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=10))
  p2 = DimPlot(obj, reduction=reduction, split.by="seurat_clusters", ncol=ncols) + theme_void() + coord_fixed(ratio=1) + NoLegend()
  plot = plot_grid(p1, p2, ncol=1, rel_heights=c(0.4,0.6))
  return(plot)
}


# RNA vs SB metrics
plot_RNAvsSB <- function(obj) {

  coords <- Misc(obj, "coords")
  if ("umi_dbscan" %in% names(coords)) {
    obj$sb_umi <- coords$umi_dbscan %>% tidyr::replace_na(0)
  } else if ("umi" %in% names(coords)) {
    obj$sb_umi <- coords$umi %>% tidyr::replace_na(0)
  } else {
    stop(trouble(
      paste0("neither 'umi_dbscan' nor 'umi' is present in the object's spatial coords table (found: ",
             paste(names(coords), collapse = ", "), ")"),
      "These columns are written by the positioning step, so this object was built by an incompatible or older version of the pipeline",
      "Re-run the spatial analysis from scratch with --spatial-analysis --force",
      "If the positioning step was interrupted, re-run --spatial-count --spatial-analysis --force"
    ), call. = FALSE)
  }

  obj$clusters <- Misc(obj,"coords")$clusters %>% tidyr::replace_na(0)
  obj$placed <- !is.na(obj$x_um) & !is.na(obj$y_um)
  
  p1 <- ggplot(obj@meta.data, aes(x=log10(nCount_RNA), y=log10(sb_umi), col=placed)) + 
    geom_point(size=0.2) + theme_bw() + xlab("log10 RNA UMI") + ylab("log10 SB UMI") + ggtitle("SB UMI vs. RNA UMI") + 
    labs(color = "placed") +
    theme(legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.background = element_blank(),
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.margin=margin(0,0,0,0,"pt"),
          legend.box.margin=margin(0,0,0,0,"pt"),
          legend.spacing.y = unit(0.1,"lines"),
          legend.key.size = unit(0.5, "lines"))
  
  d = obj@meta.data %>% rowwise %>% mutate(x=min(clusters,5)) %>% ungroup
  p2 <- ggplot(d, aes(x=as.factor(x), y=log10(nCount_RNA))) + geom_violin(scale="count") + 
    scale_x_discrete(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==5, "5+", .)}) +
    xlab("DBSCAN clusters") + ylab("log10 RNA UMI") + ggtitle("RNA UMI vs. DBSCAN cluster") + theme_classic()
  
  d = obj@meta.data %>% group_by(seurat_clusters) %>% summarize(pct.placed=paste0(round(sum(placed)/n()*100,2),"%")) %>% setNames(c("cluster","placed"))
  m = ceiling(nrow(d)/2) ; d1 = d[seq_len(m), , drop=FALSE]
  # second half: guard against nrow(d)==1 (m==1), where (m+1):nrow(d) would be 2:1 = c(2,1),
  # pulling an out-of-range NA row plus a duplicate into the table
  d2 = if (m < nrow(d)) d[(m+1):nrow(d), , drop=FALSE] else d[0, , drop=FALSE]
  p3 <- plot_grid(plot.tab(d1), plot.tab(d2), ncol=2)
  
  plot = plot_grid(gdraw("RNA vs. SB metrics"),
                   plot_grid(plot_grid(p1,p2,ncol=1), p3, ncol=2, rel_widths=c(0.5,0.5)),
                   ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}
