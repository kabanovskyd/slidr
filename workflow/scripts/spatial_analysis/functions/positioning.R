# set Lib Path
lib_path <- Sys.getenv('R_LIBS', unset = NA)
if (is.na(lib_path) || lib_path == "") {
  # fallback to default
  .libPaths()[1]
} else {
  .libPaths(lib_path)
}

# load packages
# Run positioning for each cell in input dataframe (matrix.csv -> coords.csv)
suppressMessages(library(glue)) ; g=glue ; len=length
suppressMessages(library(gridExtra))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dbscan))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Matrix))
suppressMessages(library(rdist))
suppressMessages(library(furrr))
suppressMessages(library(future))
suppressMessages(library(parallel))
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
source(paste0(r_func_path, '/positioning_func.R'))


# DBSCAN hyperparameters: search over two eps values and a range of minPts
# eps is the neighborhood radius (µm); minPts is the minimum cluster density
eps.vec = c(50, 100)
minPts.vec = c(3:42) # auto-searches up to 40x26 expansions in opt_dbscan
# KDE hyperparameters: bw controls spatial smoothing; radius defines the inclusion zone around each peak
bw = 800     # bandwidth of gaussian kernel (µm)
radius = 200 # inclusion radius around density peak (µm)

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  matrix_path <- args[[1]] # (cell, x, y, umi) dataframe from load_matrix.R
  out_path <- args[[2]]    # path to write coords.csv and diagnostic plots
  ncores <- as.numeric(args[[3]])
} else {
  stop(trouble(
    "Usage: Rscript positioning.R matrix_path output_path ncores",
    "This script is normally invoked by the pipeline's spatial-analysis stage, not by hand",
    "Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
    paste("Got", length(args), "argument(s):", paste(args, collapse = " "))
  ), call. = FALSE)
}
if (!file.exists(matrix_path)) {
  stop(trouble(
    paste0("the spatial matrix file was not found: ", matrix_path),
    "matrix.csv.gz is written by load_matrix.R in the preceding step, so this means that step did not finish",
    "Check the spatial-analysis log for load_matrix.R's own error",
    "Re-run the spatial analysis from scratch with --spatial-analysis --force"
  ), call. = FALSE)
}
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }


# load matrix.csv.gz and split into a named list of per-cell data frames
# each element is named by cb_index and contains the bead rows for that cell sorted by descending UMI
df = read.table(matrix_path, header=T, sep=",")
if (!isTRUE(all.equal(names(df), c("cb_index", "x_um", "y_um", "umi")))) {
  stop(trouble(
    paste0("unexpected columns in ", matrix_path, " (found: ", paste(names(df), collapse = ", "), ")"),
    "Expected exactly cb_index, x_um, y_um and umi, in that order",
    "A file left over from an older version of the pipeline is the usual cause",
    "Re-run the spatial analysis from scratch with --spatial-analysis --force"
  ), call. = FALSE)
}
if (nrow(df) == 0) {
  stop(trouble(
    paste0("no spatial data rows in ", matrix_path, "; there is nothing to position"),
    "No spatial-barcode read matched both a called cell and a puck bead for this sample",
    "Check the `Puck ID` metadata column names the puck this sample was actually run on",
    "Check the spatial-barcode library was sequenced deeply enough -- the read counts are in the spatial-analysis log",
    "Check the `SB Index` metadata column points at the spatial library and not the gene-expression one"
  ), call. = FALSE)
}
xlims = range(df$x_um) ; xrange = max(df$x_um) - min(df$x_um)
ylims = range(df$y_um) ; yrange = max(df$y_um) - min(df$y_um)
num_pucks = round(xrange/yrange) # estimated number of pucks, used for plot grid layout
data.list = split(df, df$cb_index)
data.list %<>% map(~arrange(.,desc(umi),cb_index))
rm(df) ; invisible(gc())
log_ts(g("positioning {len(data.list)} cells"))

### Run DBSCAN #################################################################
# opt_dbscan finds the eps and minPts that maximize % cells with exactly one cluster
plan(multisession, workers=ncores)
params = opt_dbscan(data.list)
eps = params$eps[params$is.max][[1]]
minPts = params$minPts[params$is.max][[1]]
pct.placed = round(max(params$pct)*100, 2)
log_detail(g("DBSCAN eps={eps}, minPts={minPts} → {pct.placed}% placed"))
optim_plot = ggplot(params, aes(x=minPts, y=pct*100, col=as.factor(eps))) + geom_line() +
  theme_bw() + ylab("% Placed") + labs(col="eps") + ggtitle("Parameter optimization") +
  geom_vline(xintercept = minPts, color = "red", linetype = "dashed") +
  annotate(geom = 'text', label = g("eps: {eps}\nminPts: {minPts}\nplaced: {pct.placed}%"), x = minPts+1, y = 0, hjust = 0, vjust = 0, col="red") +
  theme(legend.position = c(0.95, 0.50), legend.justification = c("right", "center"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
# apply optimal parameters to all cells; borderPoints=F excludes border points from cluster assignments
data.list %<>% lapply(function(df){
  mutate(df,
         cluster = dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$umi, borderPoints=F)$cluster,
         eps = eps,
         minPts = minPts,
         pct.placed = pct.placed)
})

dbscan_coords <- create_dbscan_coords(data.list)
invisible(gc())

plot <- plot_dbscan(dbscan_coords, optim_plot)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "DBSCAN.pdf"), 7, 8)))


### Run KDE ####################################################################
# kde() returns the top-2 density peak locations and their ratio (d2/d1) for each cell
# ratio < 1/3 means the second peak is weak relative to the dominant peak — unambiguous placement
kde_coords <- map(data.list, ~kde(., bw, radius)) %>% bind_rows
if (any(is.na(kde_coords$d1)) || any(is.na(kde_coords$d2))) {
  stop(trouble(
    "the kernel density estimation produced missing peak densities for some cells",
    "Every cell must yield two density peaks for the ambiguity ratio to be computable",
    "A cell with too few spatial beads to form a density surface is the usual cause",
    "This is an internal inconsistency rather than a configuration problem; re-run with --spatial-analysis --force",
    "If it persists after a clean re-run, please report it at https://github.com/kabanovskyd/slidr/issues"
  ), call. = FALSE)
}
kde_coords %<>% mutate(ratio = d2/d1) %>% select(1:7, ratio, everything())
plot <- plot_kde(kde_coords)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "KDE.pdf"), 7, 8)))


### More plots + save output ###################################################
if (!isTRUE(all(dbscan_coords$cb_index == kde_coords$cb_index))) {
  stop(trouble(
    "the DBSCAN and KDE coordinate tables are not in the same cell order",
    "Both are derived from the same per-cell list, so they must line up row for row before being merged",
    "This is an internal inconsistency in positioning.R rather than a problem with your data",
    "Please report it at https://github.com/kabanovskyd/slidr/issues, quoting this message"
  ), call. = FALSE)
}
coords <- merge(dbscan_coords, kde_coords, by="cb_index", suffix=c("_dbscan","_kde"))
coords %<>% rename(x2_um_kde=x2_um, y2_um_kde=y2_um)
plot <- dbscan_vs_kde(coords)
suppressMessages(suppressWarnings(make.pdf(plot, file.path(out_path, "DBSCANvsKDE.pdf"), 7, 8)))

plots <- sample_bead_plots(data.list, coords)
suppressMessages(suppressWarnings(make.pdf(plots, file.path(out_path, "beadplots.pdf"), 7, 8)))


# final placement: use the DBSCAN centroid only when ratio < 1/3 (KDE confirms one dominant location)
# cells with ratio >= 1/3 are left unplaced (x_um = y_um = NA) to avoid ambiguous assignments
coords %<>% mutate(x_um = ifelse(ratio<1/3, x_um_dbscan, NA),
                  y_um = ifelse(ratio<1/3, y_um_dbscan, NA)) %>%
                  select(cb_index, x_um, y_um, everything())
write.table(coords, file.path(out_path, "coords.csv"), sep=",", row.names=F, col.names=T, quote=F)

### Position debugging (optional) ##############################################

# plot.sb <- function(subdf) {
#   subdf %<>% arrange(umi)
#   ggplot() + coord_fixed(ratio=1, xlim=xlims, ylim=ylims) + theme_void() +
#     geom_point(data=subdf, mapping=aes(x=x_um, y=y_um, col=umi), size=2, shape=20)
# }
# plot.kde <- function(subdf) {
#   if(nrow(subdf)==0) {return(gdraw("No points"))}
#   p = Nebulosa:::wkde2d(x=subdf$x_um, y=subdf$y_um, w=subdf$umi, n=200, lims=c(xlims, ylims)) %>% {transmute(reshape2::melt(as.matrix(.[[3]])), x_um=.[[1]][Var1], y_um=.[[2]][Var2], value=value)}
#   rowmax = p[which.max(p$value),]
#   ggplot(p, aes(x=x_um, y=y_um, fill=value))+geom_tile()+coord_fixed(ratio=1) +
#     annotate("path", x=rowmax$x_um+radius*cos(seq(0,2*pi,length.out=100)), y=rowmax$y_um+radius*sin(seq(0,2*pi,length.out=100)))
# }
# plot.metadata <- function(row) {
#   row %<>% select(-x_um, -y_um)
#   row1 <- select(row, 2:10)
#   d1 <- data.frame(names(row1), round(as.numeric(row1[1,]),2)) %>% setNames(c("Data","Value"))
#   row2 <- select(row, 11:22)
#   d2 <- data.frame(names(row2), round(as.numeric(row2[1,]),2)) %>% setNames(c("Data","Value"))
#   
#   plot_grid(gdraw(g("[{row$cb_index}]")),
#             plot_grid(plot.tab(d1), plot.tab(d2), ncol=2),
#             ncol=1, rel_heights=c(0.05,0.95))
# }
# debug_coords <- function(data.list, coords) {
#   library(shiny)
#   ui <- fluidPage(
#     fluidRow(
#       column(6, plotOutput("plot1", click = "plot1_click")),
#       column(6, plotOutput("plot2"))
#     ),
#     fluidRow(
#       column(6, plotOutput("plot3")),
#       column(6, plotOutput("plot4"))
#     )
#   )
#   
#   server <- function(input, output) {
#     output$plot1 <- renderPlot({ ggplot(coords, aes(x=x_um, y=y_um)) + geom_point() + coord_fixed() + theme_void() })
#     output$plot2 <- renderPlot({ plot.new() })
#     output$plot3 <- renderPlot({ plot.new() })
#     output$plot4 <- renderPlot({ plot.new() })
#     
#     observeEvent(input$plot1_click, {
#       click <- input$plot1_click
#       if(!is.null(click)) {
#         dists <- sqrt((coords$x_um - click$x)^2 + (coords$y_um - click$y)^2)
#         row <- coords[which.min(dists),]
#         cb_index <- row$cb_index %>% as.character
#         df <- data.list[[cb_index]]
#         output$plot2 <- renderPlot({ plot.sb(df) })
#         output$plot3 <- renderPlot({ plot.kde(df) })
#         output$plot4 <- renderPlot({ plot.metadata(row) })
#       }
#     })
#   }
#   print(shinyApp(ui = ui, server = server))
#   return(T)
# }
# debug_coords(data.list, coords)
