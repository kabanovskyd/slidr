# main entry point for the spatial analysis stage
# validates input paths and launches one run_positioning.R subprocess per sample

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
suppressMessages(library(jsonlite))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(Seurat))
suppressMessages(library(rlist))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(qpdf))
suppressMessages(library(qs))
suppressMessages(library(rhdf5))
suppressMessages(library(ggrastr))
suppressMessages(library(Matrix))
suppressMessages(library(dbscan))
suppressMessages(library(rdist))
suppressMessages(library(furrr))


# Suppress warnings
options(warn = -1)
# # Restore warning settings
# options(warn = getOption("warn"))
set.seed(42)
options(future.globals.maxSize = Inf)


# R_FUNC points to the functions/ directory; source helper functions and locate run_positioning.R
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
run_positioning_src <- paste0(r_func_path, '/../run_positioning.R')


# args: BCL name, bracketed sample list (e.g. "[Sample1, Sample2]"), cellbender flag, optional ncores
# check_args_and_get_paths resolves per-sample file paths and returns a data frame with one row per sample
args <- commandArgs(trailingOnly = TRUE)
result <- check_args_and_get_paths(args)

log_section("spatial analysis")
log_ts(sprintf("%d sample(s) to position", length(result$name)))

# run run_positioning.R as a subprocess for each sample so each has its own isolated R environment
for (n in result$name) {
  t_sample <- Sys.time()
  RNAh5path = result[result$name == n, ]$rna
  Molpath = result[result$name == n, ]$molecule
  SummaryPath = result[result$name == n, ]$summary
  SBh5path = result[result$name == n, ]$spatial
  OutPath = result[result$name == n, ]$output
  ncores = result[result$name == n, ]$ncores
  PercentUmiFilter = result[result$name == n, ]$percent_umi_filter
  # shQuote every arg: system2() does NOT shell-escape args (it pastes them unquoted into
  # /bin/sh -c), so a path derived from Sample Name/BCL metadata containing shell metacharacters
  # (e.g. `;`, `$(...)`) would otherwise execute as a command
  positioning_result <- system2("Rscript", args = shQuote(c(run_positioning_src, RNAh5path, Molpath, SummaryPath, SBh5path, OutPath, as.character(ncores), as.character(PercentUmiFilter))))
  if (positioning_result != 0) {
    stop(trouble(
      paste0("run_positioning.R failed with code ", positioning_result, " for sample '", n, "'"),
      "run_positioning.R's own error is above -- it names the real cause",
      paste0("Partial outputs for this sample were left in ", OutPath),
      "Fix the cause, then re-run this stage alone with --spatial-analysis --force",
      paste0("Samples processed before '", n, "' completed successfully and do not need redoing")
    ), call. = FALSE)
  }
  log_ts(sprintf("✓ %s — %s", n, format_duration(difftime(Sys.time(), t_sample, units = "secs"))))
}
