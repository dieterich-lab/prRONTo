# TODO
# outlier
# known modification
# targets
# vicintiy
# targets, missed modifications
# ? combine barplot and table or separated
# FIXME no outliers

library(ggplot2)
library(magrittr)

OUTLIER_COL = "#1f77b4"
INLIER_COL = "#c0c0c0"

#knownmodifcation
#target
#vicinity


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
  positional_arguments = TRUE#,
  #args = c(
  #         "--feature=M",
  #         "--output=tmp",
  #         "output/results/analysis/jacusa2/preprocessed/lof/neighbors~20_contamination~0.002/cond1_vs_cond2.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$feature))
stopifnot(!is.null(opts$args))

feature <- opts$options$feature
feature_col <- paste0("feature_", feature)

result <- read.table(opts$args, sep = "\t", header = TRUE) %>%
  dplyr::mutate(is_mod = dplyr::case_when(mod != "*" ~ TRUE, .default = FALSE),
                is_outlier = dplyr::case_when(!!rlang::ensym(feature_col) == 1 ~ TRUE, .default = FALSE),
                label = dplyr::case_when(is_mod & is_outlier ~ paste(pos, mod),
                                         is_outlier ~ as.character(pos))) %>%
  split(.$seqnames)



dir.create(opts$options$output, showWarnings = FALSE)
for (seq_id in names(result)) {
  data <- result[[seq_id]]
  outlier <- data %>%
    dplyr::filter(is_outlier == TRUE)

  p_barplot <- data %>%
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
      theme(legend.position = "bottom")
  tbl <- data %>%
    head() %>%
    dplyr::mutate(coords = paste0(seqnames, ":", pos)) %>%
    dplyr::select(coords, mod, feature_col) %>%
    ggpubr::ggtexttable(rows = c(),
                        cols = c("", tail(colnames(.), n = -1)))
  p <- ggpubr::ggarrange(p_barplot, tbl,
                         labels = c("A", "B"),
                         nrow = 2)
  output <- paste0(opts$options$output, "/", seq_id, ".", opts$options$device)
  ggsave(output, p, width = 20, height = 10)
  saveRDS(p, paste0(opts$options$output, "/", seq_id, ".rds"))
  print(p)
}
