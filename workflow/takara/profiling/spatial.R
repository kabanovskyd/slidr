rm(list=ls())
############### PACKAGES ###############
library(dbscan)
library(dplyr)
library(parallel)
library(ggplot2)
library(ggthemes)
library(scales)
library(ggpubr)
library(purrr)
library(magrittr)


################ FUNCTIONS ###############
gen_NucleiPlot <- function(nuclei_toplot, minPts, plotdir, coord_table, cb_sb_umi_table){
  lapply(nuclei_toplot, function(x){
    cell_bc <- x;#print(cell_bc)
    coord_x <- coord_table[coord_table[,"cell_bc"]==cell_bc, "x_um"]
    coord_y <- coord_table[coord_table[,"cell_bc"]==cell_bc, "y_um"]
    cell <- cb_sb_umi_table %>% filter(CB == cell_bc)
    db <- dbscan::dbscan(cell[3:4], eps = eps, minPts = minPts, weights = cell$nUMI)
    cell$cluster <- db$cluster
    number_clusters <- max(cell$cluster)
    cell <- cell[order(cell$nUMI),]
    low = min(cell$nUMI);high = max(cell$nUMI)
    midpt = 0.5*(low+high)
    if(low==high){
      bounds=c(low)
    }else{
      bounds=seq(low, high, length.out=5)
    }
    centroids <- aggregate(. ~ cluster, data = cell[, c("x_um","y_um","cluster")], FUN = median)
    centroids <- centroids[centroids$cluster>0,]

    p <- ggplot(cell, aes(x = x_um, y = y_um, colour = nUMI))+
     geom_point(alpha = 0.4, size = 1)+
     theme_few()+
     xlab("")+ylab("")+
     ggtitle(paste0(cell_bc, " (", number_clusters,")"))+
     coord_fixed(ratio = 1)+
     scale_color_gradient2(low = "gray", mid = "darkblue", high = "red", midpoint = midpt,
                          limits = c(low, high),
                          breaks = bounds,
                          labels = sprintf("%.1f", bounds))+
     annotate("text", label = c("*", rep("-", nrow(centroids)))[seq_len(nrow(centroids))], x=centroids$x_um+300, y=centroids$y_um, size = 6, colour = "black", fontface = "bold")
     cell$order <- seq(nrow(cell),1, by=-1)

    q <- ggplot(cell, aes(x=order, y=nUMI))+
     geom_line()+
     geom_point()+
     theme_few()+
     xlab("Beadbarcode sorted by decreasing UMI")+
     ylab("nUMI")+
     scale_x_continuous(trans='log2', breaks = trans_breaks("log2", function(x) 2^x),
                        labels = trans_format("log2", math_format(2^.x)))+
     theme(plot.margin = margin(0.2, 0.2, 0.3, 0.3, "in"), axis.title.x=element_text(size=10))

    plt <- ggarrange(p, q, nrow = 1, ncol = 2, widths=c(3,2))
    ggsave(file.path(plotdir, paste0(x, ".jpeg")), plot = plt, width = 10, height = 6)
  })
}


generate_bcs_x_nlocations <- function(coords, n, col_name="number_clusters_o") {
  # Function to generate a list of cell barcodes based on the number of spatial locations identified
  # Args:
  #   coords: A data frame containing a column named `cell_bc` and a column named by 'col_name'
  #   n: The total number of clusters to iterate through
  # Returns:
  #   A named list of cell barcodes categorized by the number of spatial locations identified

  bcs_x_nlocations <- lapply(2:n, function(x) {
    if (x < 4) {
      bcs <- coords$cell_bc[coords[[col_name]] == x]
      print(length(bcs))
    } else {
      bcs <- coords$cell_bc[coords[[col_name]] >= x]
      print(length(bcs))
    }
    return(bcs)
  });names(bcs_x_nlocations) <- as.character(2:n)
  return(bcs_x_nlocations)
}


count_SalvagableNuclei <- function(cell_bc_list, eps, minPts, df_split){
  message("Total cell groups: ", length(cell_bc_list))

  list_sal <- lapply(cell_bc_list, function(bc_list){
    message("Cell count in this group: ", length(bc_list))

    gen_empty_df <- function() {
        data.frame(cell_bc = character(0),x_um = numeric(0),y_um = numeric(0),number_clusters = numeric(0))
    }

    count_beads <- function(cell) {
        db <- dbscan::dbscan(cell[c("x_um", "y_um")], eps = eps, minPts = minPts, weights = cell$nUMI)
        cell$cluster <- db$cluster
        nclusters <- max(db$cluster)
        cell <- cell[cell$cluster > 0, ]

        a <- sapply(seq_len(nclusters), function(i) {
          nrow(unique(cell[cell$cluster == i, c("x_um", "y_um")]))
        })
        nbeads <- sum(a == 1)

        if(nclusters > 1 && nbeads == nclusters - 1){
          # salvagable
          idx_n <- which(a!=1)
          cell_top_cluster <- cell %>% dplyr::filter(cluster == idx_n)
          xc = matrixStats::weightedMedian(cell_top_cluster$x_um,w=cell_top_cluster$nUMI)
          yc = matrixStats::weightedMedian(cell_top_cluster$y_um,w=cell_top_cluster$nUMI)
          return(c(xc, yc, 1))
        }else{
          # un-salvagable
          return(c(0,0,0))
        }
    }

    if(length(bc_list)==0){
        return(gen_empty_df())
    }else{
        df_split_x=df_split[bc_list]
        message("Cells in this group: ", length(df_split_x))

        n_salvagable <- lapply(df_split_x, count_beads) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        setNames(c("x_um", "y_um", "number_clusters")) %>%
        tibble::rownames_to_column(var = "cell_bc")

        n_salvagable$x_um <- as.numeric(n_salvagable$x_um)
        n_salvagable$y_um <- as.numeric(n_salvagable$y_um)
        n_salvagable$number_clusters <- as.numeric(n_salvagable$number_clusters)
        #print(n_salvagable)
        return(n_salvagable)
   }
  })
  names(list_sal)<-names(cell_bc_list)
  return(list_sal)
}

generate_salvage_summary <- function(bcs_x_nlocs, salvage_out, n, total_nuclei) {
  # Function to generate a summary table of nuclei, and nuclei containing beads
  # Args:
  #   bcs_x_nlocs: A list of barcodes for nuclei locations
  #   salvage_out: A list of data containing the `saved` field for each group
  #   n: The total number of clusters
  # Returns:
  #   A data frame summarizing nuclei, salvaged nuclei

  # Count nuclei based on # spatial locations
  Nuclei <- sapply(bcs_x_nlocs[as.character(2:n)], function(x) length(x))
  print(Nuclei)

  # Assign names
  if (n >= 4) {
    names(Nuclei) <- c("2", "3", ">=4")
  } else {
    names(Nuclei) <- as.character(2:n)
  }

  # Count nuclei that can be salvaged due to bead contamination
  Nuclei_with_beads <- sapply(salvage_out[as.character(2:n)], function(x) sum(x$number_clusters))
  print(Nuclei_with_beads)

  # Assign names
  if (n >= 4) {
    names(Nuclei_with_beads) <- c("2", "3", ">=4")
  } else {
    names(Nuclei_with_beads) <- as.character(2:n)
  }

  # Combine into a summary table
  sal_summary <- rbind(Nuclei, Nuclei_with_beads)

  # Ensure all expected columns ("2", "3", ">=4") are present
  expected_cols <- c("2", "3", ">=4")
  missing_cols <- setdiff(expected_cols, colnames(sal_summary))

  if (length(missing_cols) > 0) {
          for (col in missing_cols) {
                sal_summary <- cbind(sal_summary, setNames(as.data.frame(matrix(0, nrow = nrow(sal_summary), ncol = 1)), col))
          }
          # Reorder columns to expected order
          sal_summary <- sal_summary[, expected_cols, drop = FALSE]
  }
  return(sal_summary)
}


run_dbscan_sweep <- function(data_list, eps.vec, minPts_vec, cores) {
  params = expand.grid(eps.vec, minPts_vec) ; colnames(params) = c("eps", "minPts")

  start.time <- Sys.time(); print(start.time)
  mat = mclapply(1:nrow(params), function(i) {
    eps = params$eps[[i]]; minPts = params$minPts[[i]]
    m = data_list %>% map_lgl(~max(dbscan::dbscan(.[c("x_um","y_um")], eps=eps, minPts=minPts, weights=.$nUMI)$cluster) == 1)
    p = data_list %>% map_lgl(~max(dbscan::dbscan(.[c("x_um","y_um")], eps=eps, minPts=minPts, weights=.$nUMI)$cluster) >= 1)
    return(round(c(eps, minPts, 100 * (sum(m) / length(m)), 100 * (sum(p) / length(m))), 2))
  }, mc.cores = cores) %>% {do.call(rbind, .)} %>% as.data.frame ; colnames(mat) = c("eps", "minPts", "pctConfPostioned", "pctPostioned")
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time, 2); print(time.taken)
  return(mat)
}


get_minPts <- function(mat) {
  minPts = mat$minPts[which.max(mat$pctConfPostioned)]
  return(minPts)
}


rerun_dbscan <- function(mat, minPts, data_list_split, col_names) {
  start.time <- Sys.time(); print(start.time)
  #eps = mat$eps[mat$is.max][[1]]
  eps = mat$eps[which.max(mat$pctConfPostioned)]
  print(paste0("Optimal minPts is ", minPts))
  print(paste0("eps is ", eps))
  coords = lapply(data_list_split, function(df) {
    df$cluster <- dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$nUMI)$cluster
    SB_total <- length(df$matched_beadbarcode)
    SB_UMI_total <- sum(df$nUMI)
    cell_noise <- df %>% filter(cluster == 0)
    SB_noise <- length(cell_noise$matched_beadbarcode)
    SB_UMI_noise <- sum(cell_noise$nUMI)
    cell_top_cluster <- df %>% filter(cluster == 1)
    SB_top_cluster <- length(cell_top_cluster$matched_beadbarcode)
    SB_UMI_top_cluster <- sum(cell_top_cluster$nUMI)
    m = max(df$cluster) ; umimax = max(df$nUMI)
    nc_0 = m
    dmax = 0
    if (m >= 1) {
	# remove 0 cluster, and cluster with only 1 unique spatial coords (beads)
	df1 <- df %>% 
		filter(cluster > 0) %>%
		group_by(cluster) %>%
		filter(n_distinct(x_um, y_um) > 1,
			max(nUMI, na.rm = TRUE) > UMICUTOFF
		) %>%
		ungroup()
	if (nrow(df1) == 0){
		m=0 #bead
		return(c(0, 0, m, SB_total, SB_UMI_total, SB_noise, SB_UMI_noise, SB_top_cluster, SB_UMI_top_cluster, umimax, nc_0, dmax))
	} else {
		if (m>1){
			dmax <- calculate_distance(df1)$max_dist
                        	if (dmax <= DCUTOFF){
                          	m=1
                        }
                        i <- which.max(replace(df1$nUMI, is.na(df1$nUMI), -Inf))
                        cmax <- df1$cluster[i]
			df1 %<>% dplyr::filter(cluster==cmax)
		}
		xmu = matrixStats::weightedMedian(df1$x_um,w=df1$nUMI)
      		ymu = matrixStats::weightedMedian(df1$y_um,w=df1$nUMI)
	}
        return(c(xmu, ymu, m, SB_total, SB_UMI_total, SB_noise, SB_UMI_noise, SB_top_cluster, SB_UMI_top_cluster, umimax, nc_0, dmax))
    } else {
      return(c(0, 0, m, SB_total, SB_UMI_total, SB_noise, SB_UMI_noise, SB_top_cluster, SB_UMI_top_cluster, umimax, nc_0, dmax))
    }
  }) %>% do.call(rbind,.) %>% as.data.frame() %>% setNames(col_names) %>%
    tibble::rownames_to_column(var = "cell_bc")
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time, 2); print(time.taken)
  return(coords)
}

calculate_distance <- function(cell){
  centroids <- cell %>%
    filter(!is.na(cluster), nUMI > 0) %>%   # guardrails
    group_by(cluster) %>%
    summarise(
      x_um = matrixStats::weightedMedian(x_um, w = nUMI, na.rm = TRUE),
      y_um = matrixStats::weightedMedian(y_um, w = nUMI, na.rm = TRUE),
      .groups = "drop"
  )
  dist_mat <- dist(centroids[, c("x_um", "y_um")])
  max_dist <- if (length(dist_mat) == 0) 0 else max(dist_mat) 
  list(max_dist = max_dist, all_dist = dist_mat)
}

##################################################
message("########### INPUT CHECKS ###############")
##################################################
args = commandArgs(trailingOnly=TRUE);print(args);print(length(args))
if (length(args)!=6) {
  stop(paste0("expected 6 arguments: Sample ID, Input dir, Output dir, sc Platform, Subsample, Cores",
       "\nTroubleshooting:",
       "\n \u2022 This script is normally invoked by the Flex spatial-analysis stage, not by hand",
       "\n \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
       "\n \u2022 Received ", length(args), " argument(s): ", paste(args, collapse=" ")), call.=FALSE)
} else  {
  sample_id = args[1]
  in_dir = args[2] #misc
  out_dir = args[3] #output
  scPlatform = args[4]
  subsample = args[5]
  cores = as.integer(args[6])
}


##################################################
message("############ OUTPUT PATHS #############")
##################################################
coords_output_path <- file.path(out_dir, paste0("coords_", sample_id, ".txt"))
general_plots_path <- file.path(out_dir, "plots")

dir.create(general_plots_path)

##################################################
message("############## PARAMETERS #############")
##################################################
nUMI = 256
S = 0.647 # pixel to um  (10 x 10 mm trekker tile)
eps = 50

if (scPlatform %in% c("tenx","TrekkerU_C","TrekkerC","TrekkerU_CX","TrekkerCX")){
  print(scPlatform)
  minPts_vec = c(4:15)
  DCUTOFF = 0
  UMICUTOFF = 1
}else if (scPlatform %in% c("Rhapsody","TrekkerR","TrekkerU_R","TrekkerU_RATAC","TrekkerU_RVDJ")){
  print(scPlatform)
  minPts_vec = seq(4,31,by=2)
  DCUTOFF = 0
  UMICUTOFF = 1
}else if (scPlatform %in% c("TrekkerU_M")){
  print(scPlatform)
  minPts_vec = seq(4,26,by=2)
  DCUTOFF = 0
  UMICUTOFF = 1
}else if (scPlatform %in% c("TrekkerFX_FLEX")){
  print(scPlatform)
  minPts_vec = seq(4,80,by=2)
  DCUTOFF = 1000
  UMICUTOFF = 0
}else if (scPlatform %in% c("TrekkerQ_S","TrekkerQ_P","TrekkerU_IL","TrekkerU_PIP","Trekker5C_CX","Trekker5C_C")){
  print(scPlatform)
  minPts_vec = seq(4,60,by=2)
  DCUTOFF = 0
  UMICUTOFF = 0
}else if (scPlatform %in% c("TrekkerSHA_WGA")){
  print(scPlatform)
  minPts_vec = seq(4,60,by=2)
  DCUTOFF = 0
  UMICUTOFF = 0
}else {
  stop(paste0("Unrecognized single cell platform: ", scPlatform,
  ".\nPlease choose among: TrekkerU_C, TrekkerU_CX, TrekkerC, TrekkerCX, TrekkerFX_FLEX, TrekkerU_IL, TrekkerU_PIP, TrekkerSHA_WGA, TrekkerU_R, TrekkerU_RATAC, TrekkerU_RVDJ, TrekkerR, TrekkerU_M, TrekkerQ_S, TrekkerQ_P, Trekker5C_CX, Trekker5C_C",
       "\nTroubleshooting:",
       "\n \u2022 The pipeline sets this to TrekkerFX_FLEX for Flex runs, so an unexpected value means the Trekker samplesheet was hand-edited",
       "\n \u2022 Re-run the stage through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis --force`"))
}

print("range of minPts to perform sweep:") 
print(minPts_vec)


##################################################
message("############# LOAD DATA ###############")
##################################################
start.time <- Sys.time();print(start.time)
sb_whitelist <- read.delim(file.path(in_dir, paste0("df_whitelist_", sample_id, ".txt")), sep = ",")
sb_matched <- read.csv(file.path(in_dir, paste0("matching_result_", sample_id, ".csv")))
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2);print(time.taken)


##################################################
message("############# DATA SETUP ###############")
##################################################
start.time <- Sys.time();print(start.time)
# scale coordinates to microns
sb_matched$x_um <- sb_matched$xcoord*S
sb_matched$y_um <- sb_matched$ycoord*S
# clean up data frame (get rid of unscaled coordinates)
sb_matched <- sb_matched %>% dplyr::select(matched_beadbarcode, Illumina_barcode, x_um, y_um)

# merge matching results with whitelist
cb_sb_spatial_data <- merge(sb_matched, sb_whitelist, by.x = "Illumina_barcode", by.y = "SB")
# remove non-distinct rows (if there are any, there technically shouldn't be)
cb_sb_spatial_data <- cb_sb_spatial_data %>% distinct()
# cell barcode - spatial barcode combinations 
cb_sb_spatial_data <- cb_sb_spatial_data[order(cb_sb_spatial_data$nUMI), ]
## collapse illumunia_barcode into matched_beadbarcode and recalculate nUMI
cb_sb_spatial_data_collapse <- aggregate(cb_sb_spatial_data$nUMI, list(cb_sb_spatial_data$matched_beadbarcode), FUN=sum) 
names(cb_sb_spatial_data_collapse) <- c("matched_beadbarcode", "nUMI_collapsed")
coordinates <- cb_sb_spatial_data %>% dplyr::select(matched_beadbarcode, x_um, y_um)
cb_sb_spatial_data_collapse <- merge(cb_sb_spatial_data_collapse, coordinates, by = "matched_beadbarcode", all = FALSE)
cb_sb_spatial_data_collapse <- cb_sb_spatial_data_collapse[!duplicated(cb_sb_spatial_data_collapse),]


##################################################
message("#### REMOVE SB WITH VERY HIGH UMIS ####")
##################################################
cb_sb_spatial_data_collapse_filtered <- cb_sb_spatial_data_collapse %>% dplyr::filter(nUMI_collapsed < nUMI)
cb_sb_spatial_data_filtered_by_total <- merge(cb_sb_spatial_data_collapse_filtered, cb_sb_spatial_data, by = "matched_beadbarcode")
cb_sb_spatial_data_filtered_by_total <- cb_sb_spatial_data_filtered_by_total %>% dplyr::select(matched_beadbarcode,
                                                                                               nUMI_collapsed,
                                                                                               x_um.x,
                                                                                               y_um.x,
                                                                                               Illumina_barcode,
                                                                                               nUMI,
                                                                                               CB
)
names(cb_sb_spatial_data_filtered_by_total) <- c("matched_beadbarcode","nUMI_collapsed","x_um","y_um",
                                                 "Illumina_barcode","nUMI","CB")
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2);print(time.taken)

#saveRDS(cb_sb_spatial_data_filtered_by_total, paste0(sample_id, "_cb_sb_spatial_data_filtered_by_total.rds"))

##################################################
message("### RUN DBSCAN with minPts = ",minPts_vec," ###")
##################################################
print(paste0("cores=", cores))
normalize = F
rpm = F
data.list_split = split(cb_sb_spatial_data_filtered_by_total, cb_sb_spatial_data_filtered_by_total$CB)
if (normalize) {data.list_split%<>%resample.normalize}
total_cells <- length(data.list_split);print(total_cells)

# with zero cells surviving upstream filtering, both the subsample branch and the "no" branch below
# assign an empty data.list, which crashes deep inside rerun_dbscan()'s setNames() with an opaque
# "'names' attribute [N] must be the same length as the vector [0]" error -- fail clearly instead
if (total_cells == 0) {
   stop(paste0("no cells present for sample '", sample_id, "' after total-UMI filtering; cannot run spatial positioning",
       "\nTroubleshooting:",
       "\n \u2022 Every nucleus was filtered out before positioning, so there is nothing left to place",
       "\n \u2022 Check the `Puck ID` metadata column names the puck this sample was actually run on -- the wrong puck places no nuclei",
       "\n \u2022 Check the spatial-barcode library was sequenced deeply enough; the read and UMI counts are in the run's takara_pipeline.log",
       "\n \u2022 Check the `Flex Probe Barcode IDs` metadata column matches the barcodes actually used, so reads are not all discarded"), call.=FALSE)
}

set.seed(123)
if (subsample == "yes") {
   subsample_nuclei <- min(5000, floor(total_cells / 5))
   if (subsample_nuclei == 0) {
      # too few cells to take a non-empty 1/5th subsample (total_cells < 5) -- an empty
      # data.list here would silently produce NaN/numeric(0) downstream in run_dbscan_sweep()/
      # get_minPts() and crash deep inside rerun_dbscan()'s dbscan() call with an opaque
      # "argument of length zero" error, so fall back to using all cells instead
      message("Too few cells (", total_cells, ") to subsample; using all cells instead")
      data.list <- data.list_split
   } else {
      print(subsample_nuclei)
      sweep_cell_bcs <- sample(seq_along(data.list_split), size = subsample_nuclei)
      data.list <- data.list_split[sweep_cell_bcs]
      print(paste0("Nuclei Subsampled to ", subsample_nuclei))
   }
} else {
   data.list <- data.list_split
}
mat <- run_dbscan_sweep(data.list, c(50), minPts_vec, cores)
write.csv(mat, file.path(general_plots_path, "summary_minPts.csv"), row.names = FALSE, quote = FALSE)

minPts_sweep_1 <- get_minPts(mat)


####################################################
message("#### RERUN DBSCAN with optimal MINPTS ####")
####################################################
coords <- rerun_dbscan(mat, minPts_sweep_1, data.list_split, c("x_um_o","y_um_o","number_clusters_o","SB_total","SB_UMI_total","SB_noise","SB_UMI_noise","SB_top_cluster","SB_UMI_top_cluster","umimax","nc_o","dist_max"))


##################################################
message("#### CLEANUP COORDINATES DATAFRAME ####")
##################################################
start.time <- Sys.time();print(start.time)
for (col in seq_len(ncol(coords))[-1]){
  coords[,col] <- as.numeric(coords[,col])
}

## create column with proportion SB noise and proportion SB UMI noise
coords <- coords %>%
  mutate(proportion_SB_noise = SB_noise / SB_total) %>%
  mutate(proportion_SB_UMI_noise = SB_UMI_noise / SB_UMI_total) %>%
  mutate(proportion_SB_top_cluster = SB_top_cluster / SB_total) %>%
  mutate(proportion_SB_UMI_top_cluster = SB_UMI_top_cluster / SB_UMI_total) %>%
  mutate(rep_status = 0)

end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2);print(time.taken)


##################################################
message("############# BEAD REMOVAL #############")
##################################################
start.time <- Sys.time();print(start.time)

## categorize bcs into categories
#n<-min(4, max(coords$number_clusters_o))
n<-4
message("grouping nuclei by the number of spatial location assigned up to: ", n)
bcs_x_nlocs <- generate_bcs_x_nlocations(coords, n)


## count salvagable nuclei
message("eps: ", eps)
message("minPts: ", minPts_sweep_1)

salvage_out <- count_SalvagableNuclei(bcs_x_nlocs[as.character(2:n)],
                                      eps=eps,
                                      minPts=minPts_sweep_1,
                                      df_split=data.list_split)

end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2);print(time.taken)


###################################################
message("### SUMMARIZE THE BEAD REMOVAL RESULT ###")
###################################################
summary_table <- generate_salvage_summary(bcs_x_nlocs, salvage_out, n, total_nuclei=total_cells)

write.csv(summary_table, file.path(in_dir, paste0(sample_id,"_beadRemoval_summary.csv")))


###################################################
message("##### OUTPUT CELL BARCODES SALVAGED #####")
###################################################
coords_salvaged<-do.call("rbind", salvage_out)
coords_salvaged<-coords_salvaged[coords_salvaged$number_clusters==1,]
row.names(coords_salvaged)<-NULL
dim(coords_salvaged)

write.csv(coords_salvaged, file.path(in_dir, paste0(sample_id,"_beadRemoved_nuclei.csv")), row.names = FALSE, quote = FALSE)


########################################################
message("##### Update Coords with salvaged nuclei #####")
coords <- coords %>%
  left_join(coords_salvaged, by = "cell_bc")
for (col in colnames(coords)) {
  if (col %in% colnames(coords) && paste0(col, "_o") %in% colnames(coords)) {
    coords[[col]][is.na(coords[[col]])] <- coords[[paste0(col, "_o")]][is.na(coords[[col]])]
  }
}
idx <- coords$number_clusters_o != coords$number_clusters
coords$rep_status[idx] <- 1

print(scPlatform)
if (scPlatform == "TrekkerFX_FLEX") {
  ##################################################
  message("############# REP ZERO #############")
  ##################################################
  
  # get cell barcodes with 0 spatial location assignment
  cells_0 <- coords$cell_bc[coords$nc_o == 0]
  message("Number of cells with 0 spatial location assigned ", length(cells_0))
  print(length(cells_0))
  print(head(cells_0))

  if (length(cells_0) == 0) {
    message("No cells with 0 spatial locations found; skipping REP ZERO rescue pass")
  } else {
  data.list_split_zeros <- data.list_split[cells_0]

  mat_zeros <- run_dbscan_sweep(data.list_split_zeros, c(50), minPts_vec, cores)
  write.csv(mat_zeros, file.path(general_plots_path, "zeros_summary_minPts.csv"), row.names = FALSE, quote = FALSE)

  minPts_sweep_2 <- get_minPts(mat_zeros)


  ####################################################
  message("#### RERUN DBSCAN with optimal MINPTS ####")
  ####################################################
  coords_zero <- rerun_dbscan(mat_zeros, minPts_sweep_2, data.list_split_zeros, c("x_um_2","y_um_2","number_clusters_2","SB_total_2","SB_UMI_total_2","SB_noise_2","SB_UMI_noise_2","SB_top_cluster_2","SB_UMI_top_cluster_2","umimax_2","nc_o_2","dist_max_2"))


  ##################################################
  message("#### CLEANUP COORDINATES DATAFRAME ####")
  ##################################################
  start.time <- Sys.time();print(start.time)
  for (col in seq_len(ncol(coords_zero))[-1]){
    coords_zero[,col] <- as.numeric(coords_zero[,col])
  }

  ## create column with proportion SB noise and proportion SB UMI noise
  coords_zero <- coords_zero %>%
    mutate(proportion_SB_noise_2 = SB_noise_2 / SB_total_2) %>%
    mutate(proportion_SB_UMI_noise_2 = SB_UMI_noise_2 / SB_UMI_total_2) %>%
    mutate(proportion_SB_top_cluster_2 = SB_top_cluster_2 / SB_total_2) %>%
    mutate(proportion_SB_UMI_top_cluster_2 = SB_UMI_top_cluster_2 / SB_UMI_total_2)
  end.time <- Sys.time()


  ##########################################################
  message("##### Update Coords with 0 rep nuclei #####")
  ##########################################################
  coords_zero_subset <- coords_zero %>%
    filter(number_clusters_2 > 0)
  coords <- coords %>%
    left_join(coords_zero_subset, by = "cell_bc")
  for (col in colnames(coords)) {
    if (endsWith(col, "_2")) {
      original_col <- sub("_2$", "", col)
      if (original_col %in% colnames(coords)) {
        update_idx <- !is.na(coords[[col]])
        coords[[original_col]][update_idx] <- coords[[col]][update_idx]
	coords$rep_status[update_idx] <- 2
      }
    }
  }
  coords <- coords[, !endsWith(names(coords), "_2")]
  }
}

  
##################################################
message("##### WRITE COORDINATES DATAFRAME #####")
################################################## 
write.table(coords, coords_output_path)


##################################################
message("###### SUMMARIZE REP RESULTS ######")
##################################################
rep_summary_table <- table(coords$number_clusters, coords$rep_status)
rep_summary_matrix <- as.data.frame.matrix(rep_summary_table)
head(rep_summary_matrix)
for (col in c("0", "1", "2")) {
  if (!col %in% colnames(rep_summary_matrix)) {
    rep_summary_matrix[[col]] <- 0
  }
}
write.csv(rep_summary_matrix, file.path(in_dir, paste0(sample_id,"_rep_summary.csv")))


##################################################
message("###### CALCULATE QUALITY METRICS ######")
##################################################
start.time <- Sys.time();print(start.time)
## number cells mapped
spatially_mapped_top_cluster_o <- sum(coords$number_clusters_o == 1)
proportion_cells_mapped_o <- (spatially_mapped_top_cluster_o/total_cells)[1]
#spatially_mapped_any_cluster_number <- sum(coords$number_clusters_o > 0)
spatially_mapped_any_cluster_number <- sum(coords$number_clusters > 0)
spatially_mapped_top_cluster <- sum(coords$number_clusters == 1)
proportion_cells_mapped <- (spatially_mapped_top_cluster/total_cells)[1]

## sb per cell
median_unique_SB_per_cell <- median(coords$SB_total)
mean_unique_SB_per_cell <- mean(coords$SB_total)
median_SB_UMI_per_cell <- median(coords$SB_UMI_total)
mean_SB_UMI_per_cell <- mean(coords$SB_UMI_total)

## signal to noise ratios total
median_proportion_unique_SB_top_cluster_total <- median(coords$proportion_SB_top_cluster)
mean_proportion_unique_SB_top_cluster_total <- mean(coords$proportion_SB_top_cluster)
median_proportion_SB_UMI_top_cluster_total <- median(coords$proportion_SB_UMI_top_cluster)
mean_proportion_SB_UMI_top_cluster_total <- mean(coords$proportion_SB_UMI_top_cluster)

## signal to noise ratios for mapped cells
median_mapped_cells_proportion_unique_SB_top_cluster <- (coords %>% filter(number_clusters == 1) %>%
                                                           select(proportion_SB_top_cluster) %>%
                                                           summarise(across(everything(), median)))$proportion_SB_top_cluster
mean_mapped_cells_proportion_unique_SB_top_cluster <- (coords %>% filter(number_clusters == 1) %>%
                                                         select(proportion_SB_top_cluster) %>%
                                                         summarise(across(everything(), mean)))$proportion_SB_top_cluster
median_mapped_cells_proportion_SB_UMI_top_cluster <- (coords %>% filter(number_clusters == 1) %>%
                                                        select(proportion_SB_UMI_top_cluster) %>%
                                                        summarise(across(everything(), median)))$proportion_SB_UMI_top_cluster
mean_mapped_cells_proportion_SB_UMI_top_cluster <- (coords %>% filter(number_clusters == 1) %>%
                                                      select(proportion_SB_UMI_top_cluster) %>%
                                                      summarise(across(everything(), mean)))$proportion_SB_UMI_top_cluster

nUMI_cutoff_combined=nUMI
summary_metrics <- data.frame(eps,
                              minPts = minPts_sweep_1,
                              nUMI_cutoff_combined,
                              spatially_mapped_top_cluster_o,
                              spatially_mapped_any_cluster_number,
                              total_cells,
                              proportion_cells_mapped_o,
                              median_unique_SB_per_cell,
                              mean_unique_SB_per_cell,
                              median_SB_UMI_per_cell,
                              mean_SB_UMI_per_cell,
                              median_proportion_unique_SB_top_cluster_total,
                              mean_proportion_unique_SB_top_cluster_total,
                              median_proportion_SB_UMI_top_cluster_total,
                              mean_proportion_SB_UMI_top_cluster_total,
                              median_mapped_cells_proportion_unique_SB_top_cluster,
                              mean_mapped_cells_proportion_unique_SB_top_cluster,
                              median_mapped_cells_proportion_SB_UMI_top_cluster,
                              mean_mapped_cells_proportion_SB_UMI_top_cluster,
			      spatially_mapped_top_cluster,
			      proportion_cells_mapped
)


######################################################
#message("######### OUTPUT QUALITY METRICS ##########")
######################################################
write.csv(summary_metrics, file.path(general_plots_path,"summary_metrics.csv"))

