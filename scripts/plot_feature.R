library(ggplot2)
library(magrittr)

OUTLIER_COL = "#1f77b4"
INLIER_COL = "#c0c0c0"

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
  optparse::make_option(c("-d", "--device"),
                        type = "character",
                        default = "pdf",
                        help = "Plot type"),
  optparse::make_option(c("-f", "--feature"),
                        type = "character",
                        help = "Feature to plot"),
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
           "--feature=M",
           "--output=tmp",
           "output/results/analysis/jacusa2/preprocessed/lof/neighbors~20_contamination~0.002/cond1_vs_cond2.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$feature))
stopifnot(!is.null(opts$args))

feature <- opts$options$feature
feature_col <- paste0("feature_", feature)
col <- paste0("lof_outlier_" + feature)

result <- read.table(opts$args, sep = "\t", header = TRUE) %>%
  dplyr::mutate(is_mod = dplyr::case_when(mod != "*" ~ TRUE, .default = FALSE),
                is_outlier = dplyr::case_when(!!rlang::ensym(col) == 1 ~ TRUE, .default = FALSE),
                label = dplyr::case_when(is_mod & is_outlier ~ paste(pos, mod),
                                         is_outlier ~ as.character(pos))) %>%
  split(.$seqnames)

# TODO
# known modification
# targets
# neighbors
# outlier
# missed modifications

dir.create(opts$options$output, showWarnings = FALSE)
for (seq_id in names(result)) {
  data <- result[[seq_id]]
  outlier <- data %>%
    dplyr::filter(is_outlier == TRUE)

  p <- data %>%
    ggplot(aes(x = pos, y = !!rlang::ensym(feature_col), colour = is_outlier)) +
      geom_hline(aes(yintercept = 0)) +
      geom_segment(aes(x = pos,
                       y = !!rlang::ensym(feature_col),
                       xend = pos,
                       yend = !!rlang::ensym(feature_col) - !!rlang::ensym(feature_col))) +
      geom_point(data = outlier, size = 3) +
      geom_text(data = outlier,
                mapping = aes(label = label),
                hjust = 0, vjust = 0) +
      scale_colour_manual(values = c(INLIER_COL, OUTLIER_COL)) +
      labs(x = "Position [nt]", y = feature) +
      theme_bw() +
      theme(legend.position = "none")
  output <- paste0(opts$options$output, "/", seq_id, ".", opts$options$device)
  ggsave(output, p, width = 20, height = 10)
}
