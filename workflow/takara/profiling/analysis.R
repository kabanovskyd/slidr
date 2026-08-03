rm(list=ls())
library(Seurat)
library(Matrix)
library(sceasy)

options(future.globals.maxSize = 5000 * 1024^2)

############ FUNCTIONS #############
add_spatial<-function(seurato, matched_barcode){
  beads_ordered <- Cells(seurato)
  seurato[["SPATIAL"]] <- CreateDimReducObject(embeddings = as.matrix(matched_barcode)[beads_ordered,c(1,2)],
                                               key = "SPATIAL_", assay = DefaultAssay(seurato))
  coords = data.frame(x=-matched_barcode[beads_ordered, 2], y=matched_barcode[beads_ordered, 1],
                      row.names=beads_ordered, stringsAsFactors=FALSE)
  seurato@images$slice1 = new(Class = "SlideSeq", assay = "Spatial", key = "slice1_", coordinates = coords)
  return(seurato)
}

gen_markers_bycluster<-function(variable_features_cluster){
  markers_bycluster<-lapply(unique(variable_features_cluster$cluster), function(x){
    genes<-variable_features_cluster[variable_features_cluster$cluster==x,]
    genes.top<-genes[genes$p_val_adj<0.05,]
    genes.top<-genes.top[order(genes.top$avg_log2FC, decreasing = T),]
  });names(markers_bycluster)<-unique(variable_features_cluster$cluster)
  return(markers_bycluster)
}

# translate cell barcode from tenx catpure oligo to polyT oligo
translate_cellbarcode <- function(whitelist, positioned_coord){
  row.names(whitelist) <- whitelist$V2
  positioned_coord$remapped <- whitelist[positioned_coord$cell_bc,"V1"]
  # a cell_bc not found in the whitelist produces NA -> "NA-1" for every such row; if >=2 rows
  # are missing, row.names<- below would crash with "duplicate row.names are not allowed."
  # Fall back to the original (unmatched) cell_bc, which is unique per row, instead. Tag it with
  # a suffix that can never appear in a real whitelist$V1 translation, so this fallback is
  # guaranteed disjoint from every successfully-translated barcode (not just unlikely to collide).
  missing <- is.na(positioned_coord$remapped)
  positioned_coord$remapped[missing] <- paste0(positioned_coord$cell_bc[missing], "_UNMAPPED")
  positioned_coord$remapped_modified <- paste0(positioned_coord$remapped,"-1")
  row.names(positioned_coord) <- positioned_coord$remapped_modified
  return(positioned_coord)
}

# add spatial coord to seurat object
gen_matched_barcode <- function(positioned_coord){
  matched_barcode <- positioned_coord[,c("x_um","y_um")]
  colnames(matched_barcode) <- c("SPATIAL_1","SPATIAL_2")
  return(matched_barcode)
}

# add confidence of the positioning to seurat object
add_positioning_confidence <- function(seurato, positioned_coord){
  num_clusters <- positioned_coord[,c("number_clusters"), drop=FALSE]
  seurato <- AddMetaData(seurato, col.name = "number_clusters", num_clusters[Cells(seurato),"number_clusters"])
  return(seurato)
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

# run atac analysis
analyze_atac <- function(obj, res=0.8){
  DefaultAssay(obj) <- "ATAC"
  #obj <- NucleosomeSignal(obj)
  #obj <- TSSEnrichment(obj)
  obj <- FindTopFeatures(obj, min.cutoff = 5)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  nsv <- ncol(Loadings(obj, reduction="lsi"))
  ndim <- min(50, nsv)
  # seq_len(ndim)[-1] is "2:ndim" that degrades to an empty vector when ndim < 2, avoiding the
  # c(2,1)/c(2,0) reverse-range footgun on a degenerate LSI (drops the first LSI component as usual)
  atac_dims <- seq_len(ndim)[-1]
  obj <- RunUMAP(obj, reduction = "lsi", dims = atac_dims, reduction.name = "umap.atac", reduction.key = "atacUMAP_", assay ="ATAC")
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = atac_dims)
  obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3, resolution = res)

  # reorder cluster ID numerically
  clust_nm <- paste0("ATAC_snn_res.", res);print(paste0("Number of ATAC clusters identified: ", clust_nm))
  obj[[clust_nm]] <- factor(obj[[clust_nm, drop=TRUE]],
    levels = as.character(sort(as.numeric(levels(obj[[clust_nm, drop=TRUE]]))))
  )

  return(obj)
}

## build a joint neighbor graph using both assays
#analyze_rnaseq_atac <- function(obj, res=0.8){
#  npc <- ncol(Loadings(obj, reduction="pca"))
#  ndim <- min(30, npc)
#  nsv <- ncol(Loadings(obj, reduction="lsi"))
#  ndimsv <- min(50, nsv)
#  obj <- FindMultiModalNeighbors(
#    object = obj,
#    reduction.list = list("pca", "lsi"),
#    dims.list = list(1:ndim, 2:ndimsv),
#    modality.weight.name = c("RNA.weight","ATAC.weight"),
#    verbose = TRUE
#  )
#
#  # build a joint UMAP visualization
#  obj <- RunUMAP(
#    object = obj,
#    nn.name = "weighted.nn",
#    assay = c("RNA"),
#    verbose = TRUE
#  )
#  
#  obj <- FindNeighbors(obj, dims = 1:ndim)
#  obj <- FindClusters(obj, resolution = res)
#  
#  return(obj)
#}

# add metadata
add_meta <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  obj <- AddMetaData(obj, log(obj$nCount_RNA,base=10), col.name = c("log10_nCount_RNA"))
  obj <- AddMetaData(obj, log(obj$nFeature_RNA,base=10), col.name = c("log10_nFeature_RNA"))
  obj
}

# import expression tables from snRNAseq pipeline analysis
import_tenx_expression_table <- function(datadir){
  filtered_feature_bc_matrix <- Read10X(datadir)
  if (!is.list(filtered_feature_bc_matrix)) {
    counts <- filtered_feature_bc_matrix
  } else {
    counts <- filtered_feature_bc_matrix$`Gene Expression`
  }
  matched_sobj <- CreateSeuratObject(counts, project = sample_ID, min.cells = 0, min.features = 0)
  return(matched_sobj)
}

import_parse_scale_expression_table <- function(datadir){
  filtered_feature_bc_matrix <- Read10X(datadir)
  matched_sobj <- CreateSeuratObject(counts = filtered_feature_bc_matrix, project = sample_ID, min.cells = 0, min.features = 0)
  return(matched_sobj)
}

import_multiome_expression_table <- function(datadir){
  cat("Importing expression matrix and creating sobj using snRNAseq data...\n")
  filtered_feature_bc_matrix <- Read10X(datadir)
  matched_sobj <- CreateSeuratObject(counts = filtered_feature_bc_matrix$`Gene Expression`, project = sample_ID, min.cells = 0, min.features = 0,
  assay="RNA")

  cat("Creating ATAC assay and adding it to the sobj...\n")
  matched_sobj[["ATAC"]] <- CreateChromatinAssay(
    counts = filtered_feature_bc_matrix$Peaks,
    sep = c(":", "-"),
    min.cells = 0,
    min.features = 0
  )

  return(matched_sobj)
}

import_rhap_expression_table <- function(datadir, barcodedir, barcodeout){
  mtx <- readMM(file.path(datadir, sc_mtx))
  features <- read.table(file.path(datadir, sc_features), header = FALSE)
  barcodes_seq <- read.table(file.path(barcodedir, paste0(sample_ID, barcodeout)), header = FALSE)
  rownames(mtx) <- features$V1
  colnames(mtx) <- barcodes_seq$V1
  matched_sobj <- CreateSeuratObject(mtx)
  return(matched_sobj)
}


import_fx_expression_table <- function(datadir, barcodedir, barcodeout){
  mtx <- readMM(file.path(datadir, sc_mtx))
  features <- read.table(file.path(datadir, sc_features), header = FALSE)
  barcodes_seq <- read.table(file.path(barcodedir, paste0(sample_ID, barcodeout)), header = FALSE)
  rownames(mtx) <- features$V2
  colnames(mtx) <- barcodes_seq$V1
  matched_sobj <- CreateSeuratObject(mtx)
  return(matched_sobj)
}

import_rhap_vdj_expression_table <- function(datadir, barcodedir, barcodeout){
  mtx <- readMM(file.path(datadir, sc_mtx))
  features <- read.table(file.path(datadir, sc_features), header = FALSE)
  barcodes_id <- read.table(file.path(datadir, sc_barcodes), header = FALSE)
  barcodes_seq <- read.table(file.path(barcodedir, paste0(sample_ID, barcodeout)), header = FALSE)
  rownames(mtx) <- features$V1
  colnames(mtx) <- barcodes_seq$V1
  matched_sobj <- CreateSeuratObject(mtx)
  map_seq_id=setNames(data.frame(barcodes_seq, barcodes_id), c("seq", "id"))
  return(list(matched_sobj = matched_sobj, map_seq_id = map_seq_id))
}

import_fluent_expression_table <- function(datadir, barcodedir){
  mtx <- readMM(file.path(datadir, sc_mtx))
  features <- read.table(file.path(datadir, sc_features), header = FALSE)
  barcodes_seq <- read.table(file.path(datadir, sc_barcodes), header = FALSE)
  rownames(mtx) <- features$V2
  colnames(mtx) <- barcodes_seq$V1
  matched_sobj <- CreateSeuratObject(mtx)
  return(matched_sobj)
}

## TrekkerU_RATAC
import_rhap_atac_expression_table <- function(datadir, multidir, barcodedir){
  cat("Reading barcode renaming file...\n")
  barcodes_seq <- read.table(file.path(barcodedir, paste0(sample_ID, "_barcodes_perPublicFxn.tsv")), header = FALSE)
  
  # Import WTA matrix
  cat("Importing WTA (RNA) matrix data...\n")
  mtx <- readMM(file.path(datadir, "matrix.mtx.gz"))
  features <- read.table(file.path(datadir, "features.tsv.gz"), header = FALSE)
  barcodes <- read.table(file.path(datadir, "barcodes.tsv.gz"), header = FALSE)
  rownames(mtx) <- features$V1
  colnames(mtx) <- barcodes$V1
  
  cat("Creating Seurat object from WTA matrix...\n")
  matched_sobj <- CreateSeuratObject(mtx)
  cat("First few cell names in Seurat object:\n")
  print(head(Cells(matched_sobj)))
  
  # Import ATAC matrix
  cat("Importing ATAC matrix data...\n")
  counts <- Matrix::readMM(file.path(multidir, "atac-matrix.mtx.gz"))
  features <- read.table(file.path(multidir, "atac-features.tsv.gz"), header = FALSE)
  barcodes <- read.table(file.path(multidir, "atac-barcodes.tsv.gz"), header = FALSE)
  
  cat("First few ATAC barcodes:\n")
  print(head(barcodes))
  rownames(counts) <- features$V1
  colnames(counts) <- barcodes$V1
  
  cat("Creating ChromatinAssay for ATAC data...\n")
  atac_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    min.cells = 0,
    min.features = 0
  )
  
  cat("First few cell names in ATAC assay:\n")
  print(head(Cells(atac_assay)))
  
  cat("Adding ATAC assay to Seurat object...\n")
  matched_sobj[["ATAC"]] <- atac_assay
  
  cat("Renaming cells using barcodes sequence file...\n")
  matched_sobj <- RenameCells(matched_sobj, new.names = barcodes_seq$V1)
  
  cat("Iport complete. Returning Seurat object.\n")
  return(matched_sobj)
}


# import spatial coordinates
import_spatial_coord <- function(assay, outdir, scriptsdir){
  coord <- read.table(file.path(outdir, paste0("coords_", sample_ID,".txt")))
  if (assay %in% c("TrekkerC")){
    tenx_whitelist <- read.table(file.path(scriptsdir, "cellbarcode_whitelists", "Aug_2023_TrekkerC_3M-february-2018.txt.gz"))
    coord <- translate_cellbarcode(tenx_whitelist, coord)
    return(coord)
  }else if (assay %in% c("TrekkerCX")){
    tenx_whitelist <- read.table(file.path(scriptsdir, "cellbarcode_whitelists", "Apr_2024_TrekkerCX_3M-3pgex-may-2023.txt.gz"))
    coord <- translate_cellbarcode(tenx_whitelist, coord)
    return(coord)
  }else if (assay %in% c("TrekkerR","TrekkerU_R","TrekkerU_RATAC","TrekkerU_RVDJ","TrekkerQ_S","TrekkerQ_P","TrekkerU_IL","TrekkerU_PIP","TrekkerSHA_WGA")){
    rownames(coord) <- coord$cell_bc
    return(coord)
  }else if (assay %in% c("TrekkerU_C","TrekkerU_CX","TrekkerU_M","TrekkerFX_FLEX","Trekker5C_C","Trekker5C_CX")){
    rownames(coord) <- paste0(coord$cell_bc, "-1")
    return(coord)
  }else{
    stop(paste0("Unrecognized single cell platform: ", assay,
    ".\nPlease choose among: TrekkerU_R, TrekkerU_RATAC, TrekkerU_RVDJ, TrekkerU_C, TrekkerU_CX, TrekkerCX, TrekkerC, TrekkerFX_FLEX, TrekkerR, TrekkerU_M, TrekkerQ_S, TrekkerQ_P, TrekkerU_IL, TrekkerU_PIP, TrekkerSHA_WGA, Trekker5C_C, Trekker5C_CX",
       "\nTroubleshooting:",
       "\n \u2022 The pipeline sets this to TrekkerFX_FLEX for Flex runs, so an unexpected value means the Trekker samplesheet was hand-edited",
       "\n \u2022 Re-run the stage through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis --force`"))
  }
}

# summarize positioning rate
summarize_positioning <- function(sobj, nnuclei){
  counts_clusters <- table(sobj$number_clusters)
  with_invalid_spatial_barcodes = as.integer(nnuclei)-sum(counts_clusters)
  counts_clusters <- append(with_invalid_spatial_barcodes , counts_clusters)
  names(counts_clusters)[1] <- -1
  # consolidate positionning confidence
  if(max(as.numeric(names(counts_clusters)))>3){
    message("Max number of spatial coords identified per nuclei is at least 4.\nConsolidating nulei with at least 4...")
    counts_clusters_c <- counts_clusters[as.numeric(names(counts_clusters))<=3]
    counts_clusters_c[">=4"]<-sum(counts_clusters[as.numeric(names(counts_clusters))>3])
  }else{
    message("Max number of spatial coord identified per nuclei is fewer than 4")
    counts_clusters_c <- counts_clusters
    counts_clusters_c[">=4"]<-0
  }
  names(counts_clusters)[1] <- "with_invalid_spatial_barcodes"
  names(counts_clusters_c)[1] <- "with_invalid_spatial_barcodes"
  return(list(counts_clusters, counts_clusters_c))
}

# add vdj metadata
add_vdj_metadata <- function(seurat_obj, meta_data, vdj_columns) {
  # Check that the specified columns exist in the metadata
  missing_columns <- setdiff(vdj_columns, colnames(meta_data))
  if (length(missing_columns) > 0) {
    stop(paste0("the following columns are missing from the provided metadata: ", paste(missing_columns, collapse = ", "),
       "\nTroubleshooting:",
       "\n \u2022 These columns come from the single-cell VDJ metadata handed to this step",
       "\n \u2022 Check the upstream single-cell run produced a VDJ metadata table with all of them",
       "\n \u2022 Columns present: ", paste(colnames(meta_data), collapse = ", ")), call.=FALSE)
  }
  # Display dimensions of metadata before modification
  message("Dimensions of metadata before modification: ", paste(dim(seurat_obj@meta.data), collapse = " x "))
  # Add the VDJ metadata
  seurat_obj <- AddMetaData(seurat_obj, col.name = vdj_columns, metadata = meta_data[Cells(seurat_obj), vdj_columns])
  # Display dimensions of metadata after modification
  message("Dimensions of metadata post modification: ", paste(dim(seurat_obj@meta.data), collapse = " x "))
  return(seurat_obj)
}


#################################################################
message("################### INPUT CHECKS #####################")
#################################################################
args = commandArgs(trailingOnly=TRUE);print(args);print(length(args))
if (length(args)!=9) {
  stop(paste0("expected 9 arguments: sample ID, single cell assay platform, scripts folder, ",
       "folder containing single cell expression matrix, misc folder, output folder, ",
       "profile conda path, scmulti_outdir, sample_ID_sc",
       "\nTroubleshooting:",
       "\n \u2022 This script is normally invoked by the Flex spatial-analysis stage, not by hand",
       "\n \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
       "\n \u2022 Received ", length(args), " argument(s): ", paste(args, collapse=" ")), call.=FALSE)
} else  {
  sample_ID = args[1]
  scPlatform = args[2]
  scriptsdir = args[3]
  scdir = args[4]
  miscdir = args[5]
  outdir = args[6]
  profile_conda_path = args[7] 
  scmulti_outdir = args[8]
  sample_ID_sc = args[9]
}

source(file.path(scriptsdir, "..", "shared_seurat_utils.R"))

intermediates_dir <- file.path(outdir,"intermediates")
dir.create(intermediates_dir)

library(reticulate)
reticulate::use_condaenv(profile_conda_path, required = TRUE)

sc_barcodes <- "barcodes.tsv.gz"
sc_features <- "features.tsv.gz"
sc_mtx <- "matrix.mtx.gz"


if (scPlatform %in% c("TrekkerU_IL")) {
  sc_barcodes <- file.path("trekkerinterim", "barcodes.tsv.gz")
  sc_features <- file.path("trekkerinterim", "features.tsv.gz")
  sc_mtx <- file.path("trekkerinterim", "matrix.mtx.gz")
}

if (scPlatform %in% c("TrekkerU_M", "TrekkerU_RATAC")) {
  library(Signac)
  library(GenomicRanges)
}

#################################################################
message("### Import expression table from snRNAseq analysis ###")
#################################################################
if (scPlatform %in% c("TrekkerU_C","TrekkerC","tenx","TrekkerU_CX", "TrekkerCX" ,"Trekker5C_CX", "Trekker5C_C")){
  matched_sobj <- import_tenx_expression_table(scdir)
}else if(scPlatform %in% c("TrekkerU_M")){
  matched_sobj <- import_multiome_expression_table(scdir)
}else if (scPlatform %in% c("TrekkerU_R","TrekkerR","Rhapsody")){
  matched_sobj <- import_rhap_expression_table(scdir, miscdir, "_barcodes_perPublicFxn.tsv")
}else if (scPlatform %in% c("TrekkerFX_FLEX")){
  matched_sobj <- import_fx_expression_table(scdir, miscdir, "_barcodes_trekkerFX.tsv.gz")
}else if (scPlatform %in% c("TrekkerU_RVDJ")){
  rhap_vdj_data <- import_rhap_vdj_expression_table(scdir, miscdir, "_barcodes_perPublicFxn.tsv")
  matched_sobj <- rhap_vdj_data$matched_sobj
  map_seq_id <- rhap_vdj_data$map_seq_id
}else if(scPlatform %in% c("TrekkerU_IL", "TrekkerU_PIP", "TrekkerSHA_WGA")){
  matched_sobj <- import_fluent_expression_table(scdir)
}else if(scPlatform %in% c("TrekkerQ_S","TrekkerQ_P")){
matched_sobj <- import_parse_scale_expression_table(scdir)
}else if(scPlatform %in% c("TrekkerU_RATAC")){
  matched_sobj <- import_rhap_atac_expression_table(scdir, scmulti_outdir, miscdir)
}else{
  stop(paste0("Unrecognized single cell platform: ", scPlatform,
  ".\nPlease choose among: TrekkerU_C, TrekkerU_CX, TrekkerCX, TrekkerC, TrekkerFX_FLEX, TrekkerU_IL, TrekkerU_PIP, TrekkerSHA_WGA, TrekkerU_R, TrekkerU_RATAC, TrekkerU_RVDJ, Trekker5C_CX, Trekker5C_C, TrekkerR, TrekkerU_M, TrekkerQ_S, TrekkerQ_P",
       "\nTroubleshooting:",
       "\n \u2022 The pipeline sets this to TrekkerFX_FLEX for Flex runs, so an unexpected value means the Trekker samplesheet was hand-edited",
       "\n \u2022 Re-run the stage through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis --force`"))
}

TotalNuclei <- length(Cells(matched_sobj))
message("Total nuclei from snRNAseq pipeline analysis: ", TotalNuclei)

#################################################################
message("############## Import spatial coordinates #############")
#################################################################
coord <- import_spatial_coord(scPlatform, outdir, scriptsdir)
matched_barcode <- gen_matched_barcode(coord)
print(head(matched_barcode))

#################################################################
message("########## Combine spatial with expression ############")
#################################################################
# combine spatial with expression
matched_sobj <- subset(matched_sobj, cells = rownames(matched_barcode))
matched_sobj <- add_spatial(matched_sobj, matched_barcode)
matched_sobj <- add_positioning_confidence(matched_sobj, coord)
rm("coord")

#################################################################
message("########## Summarize positioning confidence ###########")
#################################################################
nclusters_summary <- summarize_positioning(matched_sobj, TotalNuclei)
write.table(t(as.data.frame(nclusters_summary[[1]])), file.path(miscdir, paste0(sample_ID,"_summary_position_conf.txt")), quote=FALSE, sep="\t", row.names = FALSE)
write.table(t(as.data.frame(nclusters_summary[[2]])), file.path(miscdir, paste0(sample_ID,"_summary_position_conf_consolidated.txt")), quote=FALSE, sep="\t", row.names = FALSE)

##################################################################
message("##################### Run analysis #####################")
##################################################################
matched_sobj <- add_meta(matched_sobj)
matched_sobj <- analyze(matched_sobj, res=0.2)
variable_features_cluster <- FindAllMarkers(matched_sobj, assay = "SCT", only.pos =T)

if (scPlatform == "TrekkerU_RVDJ") {
  #################################################################
  message("################### Run VDJ analysis ###################")
  message("################### VDJ Incorporation ################")
  #################################################################
  sobj_vdj <- readRDS(file.path(scmulti_outdir, paste0(sample_ID_sc, "_Seurat.rds")))
  meta_vdj=sobj_vdj@meta.data
  meta_vdj$id=rownames(meta_vdj)
  meta_vdj=merge(meta_vdj, map_seq_id, by = "id", all.x = TRUE)
  rownames(meta_vdj)=meta_vdj$seq
  col_vdj <- c("Total_VDJ_Read_Count", "Total_VDJ_Molecule_Count", "BCR_Heavy_V_gene_Dominant", "BCR_Heavy_D_gene_Dominant", "BCR_Heavy_J_gene_Dominant",
  "BCR_Heavy_C_gene_Dominant", "BCR_Heavy_CDR3_Nucleotide_Dominant", "BCR_Heavy_CDR3_Translation_Dominant", "BCR_Heavy_Read_Count", "BCR_Heavy_Molecule_Count",
  "BCR_Light_V_gene_Dominant", "BCR_Light_J_gene_Dominant", "BCR_Light_C_gene_Dominant", "BCR_Light_CDR3_Nucleotide_Dominant", "BCR_Light_CDR3_Translation_Dominant",
  "BCR_Light_Read_Count", "BCR_Light_Molecule_Count", "TCR_Alpha_Gamma_V_gene_Dominant", "TCR_Alpha_Gamma_J_gene_Dominant", "TCR_Alpha_Gamma_C_gene_Dominant",
  "TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant", "TCR_Alpha_Gamma_CDR3_Translation_Dominant", "TCR_Alpha_Gamma_Read_Count", "TCR_Alpha_Gamma_Molecule_Count",
  "TCR_Beta_Delta_V_gene_Dominant", "TCR_Beta_Delta_D_gene_Dominant", "TCR_Beta_Delta_J_gene_Dominant", "TCR_Beta_Delta_C_gene_Dominant",
  "TCR_Beta_Delta_CDR3_Nucleotide_Dominant","TCR_Beta_Delta_CDR3_Translation_Dominant", "TCR_Beta_Delta_Read_Count", "TCR_Beta_Delta_Molecule_Count",
  "BCR_Paired_Chains","TCR_Paired_Chains")
  col_vdj <- intersect(col_vdj, colnames(meta_vdj))
  matched_sobj<-add_vdj_metadata(seurat_obj=matched_sobj, meta_data=meta_vdj, vdj_columns=col_vdj)
}

if (scPlatform %in% c("TrekkerU_M", "TrekkerU_RATAC")) { 
   ##################################################################
   message("################## Run ATAC analysis ##################")
   ##################################################################
   matched_sobj <- analyze_atac(matched_sobj, res=0.8)
}

##################################################################
message("######### Subset sobj and run spatial analysis ########")
##################################################################
message("# Generate Positioned and Conf Positioned Seurat objects #")
matched_sobj_pos <- subset(matched_sobj, cells=Cells(matched_sobj)[matched_sobj$number_clusters>0])
message("Number of postioned nuclei: ", length(Cells(matched_sobj_pos)))
if (length(Cells(matched_sobj_pos)) == 0) {
  stop(paste0("no nuclei passed QC filtering (number_clusters>0) for sample ", sample_ID,
       "\nTroubleshooting:",
       "\n \u2022 No nucleus was assigned any spatial location at all",
       "\n \u2022 Check the `Puck ID` metadata column names the puck this sample was actually run on -- the wrong puck places no nuclei",
       "\n \u2022 Check the spatial-barcode library was sequenced deeply enough; the read and UMI counts are in the run's takara_pipeline.log",
       "\n \u2022 Check the `Flex Probe Barcode IDs` metadata column matches the barcodes actually used, so reads are not all discarded"), call.=FALSE)
}

matched_sobj_pos_conf <- subset(matched_sobj_pos, cells=Cells(matched_sobj_pos)[matched_sobj_pos$number_clusters==1])
message("Number of confidently positioned nuclei: ", (length(Cells(matched_sobj_pos_conf))))
if (length(Cells(matched_sobj_pos_conf)) == 0) {
  stop(paste0("no nuclei passed QC filtering (number_clusters==1) for sample ", sample_ID,
       "\nTroubleshooting:",
       "\n \u2022 Nuclei were placed, but none of them at a single unambiguous location",
       "\n \u2022 A `number_clusters==1` filter keeps only unambiguously placed nuclei, so a sample with diffuse spatial signal can pass the previous check and fail this one",
       "\n \u2022 Check the `Puck ID` metadata column names the puck this sample was actually run on -- the wrong puck places no nuclei",
       "\n \u2022 Check the spatial-barcode library was sequenced deeply enough; the read and UMI counts are in the run's takara_pipeline.log",
       "\n \u2022 Check the `Flex Probe Barcode IDs` metadata column matches the barcodes actually used, so reads are not all discarded"), call.=FALSE)
}

## re-find spatially variable features on conf positioned nuclei
DefaultAssay(matched_sobj_pos_conf)<-"SCT"
matched_sobj_pos_conf <- FindSpatiallyVariableFeatures(matched_sobj_pos_conf, assay="SCT", slot = "scale.data",
                                                      features = head(VariableFeatures(matched_sobj_pos_conf), 200),
                                                      selection.method = "moransi", x.cuts = 100, y.cuts = 100,
                                                      verbose = TRUE, nfeatures=200)
variable_features_moransi <- organize_spatial_features_moransi(matched_sobj_pos_conf)

###################################################################
message("##################### OUTPUT ###########################")
###################################################################
message("################# Save Seurat objects ##################")
DefaultAssay(matched_sobj) <- "SCT"
DefaultAssay(matched_sobj_pos) <- "SCT"
saveRDS(matched_sobj, file.path(intermediates_dir, paste0(sample_ID, "_seurat_spatial.rds")))
saveRDS(matched_sobj_pos, file.path(intermediates_dir, paste0(sample_ID, "_Positioned_seurat_spatial.rds")))
saveRDS(matched_sobj_pos_conf, file.path(outdir, paste0(sample_ID, "_ConfPositioned_seurat_spatial.rds")))

message("################ Save matched H5AD objects #############")
### adjust spatial coord before converting to h5ad object
matched_sobj <- flip_spatial(matched_sobj)
matched_sobj_pos <- flip_spatial(matched_sobj_pos)
matched_sobj_pos_conf <- flip_spatial(matched_sobj_pos_conf)
sceasy::convertFormat(matched_sobj, from="seurat", to="anndata", outFile=file.path(intermediates_dir, paste0(sample_ID, "_anndata_matched.h5ad")))
sceasy::convertFormat(matched_sobj_pos, from="seurat", to="anndata", outFile=file.path(intermediates_dir, paste0(sample_ID, "_Positioned_anndata_matched.h5ad")))
sceasy::convertFormat(matched_sobj_pos_conf, from="seurat", to="anndata", outFile=file.path(outdir, paste0(sample_ID, "_ConfPositioned_anndata_matched.h5ad")))

message("################## Save sparse matrices ################")
matched_coords <- as.data.frame(matched_sobj[["SPATIAL"]]@cell.embeddings)
output_sparse_matrix(sample_ID, matched_sobj, matched_coords, "", intermediates_dir)
output_sparse_matrix(sample_ID, matched_sobj_pos, matched_coords, "_PositionedNuclei", intermediates_dir)
output_sparse_matrix(sample_ID, matched_sobj_pos_conf, matched_coords, "_ConfPositionedNuclei", outdir)

message("################## Save variable genes #################")
write.table(variable_features_moransi, file.path(outdir, paste0(sample_ID, "_variable_features_spatial_moransi.txt")), quote = F, sep="\t")
write.csv(variable_features_cluster, file.path(outdir, paste0(sample_ID, "_variable_features_clusters.csv")), quote=F)

