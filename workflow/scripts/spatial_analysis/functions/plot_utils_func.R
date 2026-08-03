#########################################################################################
################################ SHARED PLOT HELPERS ###################################
#########################################################################################

# Format a fatal-error message followed by actionable troubleshooting bullets, in the same
# "<message> / Troubleshooting: / • ..." layout the Python side of the pipeline uses. Intended to
# wrap a stop() message: stop(trouble("what went wrong", "hint", "hint")).
#
# The spatial analysis stage runs as a chain of Rscript subprocesses whose captured output
# (spatial_analysis.log) is the only record of a failure -- the Python caller only ever sees a
# non-zero exit code -- so a stop() here has to say what to do about it, not just what went wrong.
#
# Defined in this file because it is the one helper file every other *_func.R sources.
trouble <- function(msg, ...) {
  hints <- c(...)
  if (length(hints) == 0) return(msg)
  paste0(msg, "\nTroubleshooting:\n", paste0(" • ", hints, collapse = "\n"))
}

# Progress-logging vocabulary, deliberately identical to the Python and Julia sides (helpers.log_ts /
# log_detail / log_section / format_duration, and spatial_count.jl's logts / logdetail / logsection).
# The R scripts' output is captured into spatial_analysis.log, and a reader moving between that file
# and runtime.log should not have to adjust to a different layout: one timestamp column, details
# hanging beneath it, a rule per phase.
#
# Everything goes to stderr via message(), which is where R's own diagnostics go and what the pipeline
# captures -- print()/cat() would interleave unpredictably with it.
TS_WIDTH <- 10

log_ts <- function(msg) {
  message(format(Sys.time(), "%H:%M:%S"), "  ", msg)
}

log_detail <- function(msgs) {
  for (msg in msgs) message(strrep(" ", TS_WIDTH), msg)
}

log_section <- function(title, width = 74) {
  prefix <- paste0("── ", title, " ")
  message("")
  message(prefix, strrep("─", max(0, width - nchar(prefix))))
}

# Render an elapsed time compactly: '45s', '12m34s', '2h07m' -- the same shapes the Python and Julia
# sides produce.
format_duration <- function(seconds) {
  total <- round(as.numeric(seconds))
  hours <- total %/% 3600
  minutes <- (total %% 3600) %/% 60
  secs <- total %% 60
  if (hours > 0) return(sprintf("%dh%02dm", hours, minutes))
  if (minutes > 0) return(sprintf("%dm%02ds", minutes, secs))
  sprintf("%ds", secs)
}

# Time a phase: log a start line, run `expr`, then log a tick with the elapsed time.
log_timed <- function(label, expr) {
  log_ts(label)
  started <- Sys.time()
  result <- force(expr)
  log_ts(sprintf("✓ %s — %s", label, format_duration(difftime(Sys.time(), started, units = "secs"))))
  result
}

# Helper methods
gdraw <- function(text, s=14) {ggdraw()+draw_label(text, size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df,rows=NULL)))}
add.commas <- function(num){prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
  suppressMessages(suppressWarnings({
    pdf(file = name, width = w, height = h)
    lapply(plots, function(x) { print(x) })
    garbage <- dev.off()
  }))
}
