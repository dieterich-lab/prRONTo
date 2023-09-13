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
           "output/results/analysis/jacusa2/preprocessed/lof/neighbors~20_contamination~0.002/cond1_vs_cond2.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t")
features <- colnames(df) %>%
  stringr::str_match_all("lof_outlier_(.+)$") %>%
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
  dplyr::select(seqnames, pos, strand, modified, analysis, parameters, comparison, feature_lof_outlier_cols) %>%
  tidyr::pivot_longer(cols = feature_lof_outlier_cols, names_to = "feature", values_to = "outlier") %>%
  dplyr::mutate(feature = as.factor(gsub("lof_outlier_", "", feature))) %>%
  dplyr::filter(outlier == 1) %>%
  dplyr::group_by(analysis, parameters, comparison, feature) %>%
  dplyr::summarise(outliers = dplyr::n(), mods = sum(modified)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols = c("outliers", "mods"), names_to = "count_type", values_to = "count") %>%
  dplyr::mutate(count_type = factor(count_type, c("outliers", "mods"), ordered = TRUE))

p <- filtered %>%
  ggplot(aes(x = feature, y = count, fill = count_type)) +
    geom_col(position = position_dodge()) +
    scale_fill_discrete(labels = c("outliers" = "outlier",
                                   "mods" = "outlier and modified")) +
    labs(x = "features", y = "count", fill = "") +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave(p, opts$options$output)
