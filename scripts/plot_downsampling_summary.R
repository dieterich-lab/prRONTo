library(ggplot2)
library(magrittr)

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

option_list <- list(
  optparse::make_option(c("-f", "--features"),
                        type = "character",
                        help = "','-separted list of features, otherwise guessed from column names"),
  optparse::make_option(c("-t", "--targets"),
                        type = "character",
                        help = "','-separted list of targets, e.g.: 18S:20[-22]"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE,
  args = c(
           #"--features=M",
           "--output=tmp",
           "output/results/merged_lof.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t")
features <- colnames(df) %>%
  stringr::str_match_all("lof_outlier_(.+)") %>%
  lapply(function(m) { data.frame(feature = m[, 2]) }) %>%
  unlist(use.names = FALSE)
if (!is.null(opts$options$features)) {
  features <- strsplit(opts$options$features, ",")
}
stopifnot(!is.null(features) && length(features) > 0)

# TODO
# keep targets

feature_lof_outlier_cols = paste0("lof_outlier_", features)
filtered <- df %>%
  dplyr::mutate(modified = mod != "*") %>%
  dplyr::select(seqnames, pos, strand, modified, analysis, parameters, comparison, dplyr::all_of(feature_lof_outlier_cols)) 


