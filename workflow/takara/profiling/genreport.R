#!/usr/bin/env Rscript

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

print(args)
print(length(args))


if (length(args)!=11) {
  stop(paste0("expected 11 arguments: sample_id, path to spatial barcode whitelist, htmlrender.Rmd file, ",
       "output dir, single cell assay platform, metrics file, cluster driving gene list, ",
       "spatially variable gene list, seurat obj for all nuclei, seurat object for all positioned nuclei, ",
       "mismatch frequency file",
       "\nTroubleshooting:",
       "\n \u2022 This script is normally invoked by the Flex spatial-analysis stage, not by hand",
       "\n \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
       "\n \u2022 Received ", length(args), " argument(s): ", paste(args, collapse=" ")), call.=FALSE)
} else  {

  sample = args[1] 	
  file_whitelist <- args[2] #/path/to/*_BeadBarcodes.txt
  beadbarcode_file <- strsplit(file_whitelist, split="/")[[1]][length(strsplit(file_whitelist, split="/")[[1]])]
  beadbarcode <- gsub("_BeadBarcodes.txt", "", beadbarcode_file)

  rmarkdown::render( 

	  input  = args[3], #htmlrender.Rmd,
	  output_file = if (grepl("htmlrender_extended\\.Rmd$" ,args[3])) {
		paste0(sample, "_", "Report.html")
	  } else {
		paste0(sample, "_", "Trekker_Report.html")
	  },
	  output_dir = args[4],
	  intermediates_dir = args[4],
	  knit_root_dir = args[4],	
	  params = list( 
	    sample = sample,
	    beadbarcode = beadbarcode,
	    platform = args[5],
	    file_metrics = args[6], 
	    file_variable_features_clusters = args[7],
	    file_variable_features_spatial_moransi = args[8],
	    file_matched_seurato = args[9],
	    file_matched_seurato_cpositioned = args[10],
	    file_freq_mismatch = args[11]
	)
 )
}
