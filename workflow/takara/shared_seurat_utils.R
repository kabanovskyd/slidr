library(Seurat)
library(Matrix)

# adjust coord for scanpy
flip_spatial <- function(seurato){
  coord = data.frame("SPATIAL_1"=seurato[["SPATIAL"]]@cell.embeddings[,1],"SPATIAL_2"=-seurato[["SPATIAL"]]@cell.embeddings[,2])
  beads_ordered = Cells(seurato)
  seurato[["SPATIAL"]] = CreateDimReducObject(embeddings = as.matrix(coord)[beads_ordered,c(1,2)],key = "SPATIAL_", assay = DefaultAssay(seurato))
  seurato
}

# organize spatially variable features obtained via moransi
organize_spatial_features_moransi <- function(seurato){
  variable_features_moransi <- seurato@assays$SCT@meta.features
  variable_features_moransi <- variable_features_moransi[,grep("^[Mm]orans",colnames(variable_features_moransi))]
  variable_features_moransi <- variable_features_moransi[!is.na(variable_features_moransi$MoransI_observed),]
  variable_features_moransi <- variable_features_moransi[order(variable_features_moransi$moransi.spatially.variable.rank),]
  return(variable_features_moransi)
}

# output expression tables in sparse matrix format
output_sparse_matrix <- function(sample_id, sobj, coords, suffix1, outputDir){
  spatial_coord=coords[row.names(coords)%in%Cells(sobj),]
  Matrix::writeMM(sobj[["RNA"]]@counts, file.path(outputDir, paste0(sample_id, "_", "MoleculesPer", suffix1, ".mtx")))
  write(x = rownames(sobj[["RNA"]]@counts), file = file.path(outputDir, paste0(sample_id, "_", "genes", suffix1, ".tsv")))
  write(x = colnames(sobj[["RNA"]]@counts), file = file.path(outputDir, paste0(sample_id, "_", "barcodes", suffix1, ".tsv")))
  write.csv(spatial_coord, file.path(outputDir, paste0(sample_id, "_", "Location", suffix1, ".csv")), quote = F)
}
