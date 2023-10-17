options(error = function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    sink(stderr())
    on.exit(sink(NULL))
    cat("Backtrace:\n")
    calls <- rev(calls[-length(calls)])
    for (i in seq_along(calls)) {
      cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
    }
  }
  if (!interactive()) {
    q(status = 1)
  }
})

save_plot <- function(p, fname, suffix = "", ...) {
  rds_fname <- paste0(fname, ".rds")
  if (suffix != "") {
    rds_fname <- paste0(gsub(suffix, "", fname), "rds")
  }
  saveRDS(p, rds_fname)
  ggsave(fname, p, ...)
}


read_merged_lof <- function(fname) {
  df <- read.table(opts$args, header = TRUE, sep = "\t")

  downsampling <- stringr::str_match(df$parameters, "seed~(.+)_reads~([0-9]+)") %>%
    as.data.frame()
  colnames(downsampling) <- c("parameters", "seed", "reads")
  downsampling <- downsampling %>%
    dplyr::select(-parameters)
  df <- dplyr::bind_cols(df, downsampling)

  df
}


debug_opts <- function(option_list, args) {
  if (Sys.getenv("PRONTO_DEBUG") == "TRUE") {
    write(paste0(rep("*", 80), collapse = ""), stderr())
    write("* DEBUG! IGNORING CLI, because env PRONTO_DEBUG = TRUE", stderr())
    write(paste0(rep("*", 80), collapse = ""), stderr())
    opts <- optparse::parse_args(
      optparse::OptionParser(option_list = option_list),
      positional_arguments = TRUE,
      args = args
    )
  } else {
    opts <- optparse::parse_args(
      optparse::OptionParser(option_list = option_list),
      positional_arguments = TRUE
    )
  }

  opts
}

