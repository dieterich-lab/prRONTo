# TODO
# * outlier
# * vicintiy
# * targets, missed modifications
# * FIXME no outliers
# * knownmodifcation

library(ggplot2)
library(magrittr)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

OUTLIER_COL = "#1f77b4"
INLIER_COL = "#c0c0c0"

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
  optparse::make_option(c("-n", "--neighboor"),
                        type = "integer",
                        default = 0,
                        help = "Distance to known mod up/downstream"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
args = c("--feature=M",
         "--output=tmp",
         "--targets=18S:100",
         "output/results/analysis/jacusa2/preprocessed/lof/neighbors~20_contamination~0.002/cond1_vs_cond2.tsv")
opts <- debug_opts(option_list, args)

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

# * targets
targets <- NULL
if (!is.null(opts$options$targets)) {
  targets <- strsplit(opts$options$targets, ",") %>%
    unlist() %>%
    GenomicRanges::GRanges()
}

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

  tbl = "no data"
  if (nrow(outlier) > 0) {
    tbl <- outlier %>%
      dplyr::mutate(coords = paste0(seqnames, ":", pos)) %>%
      dplyr::select(coords, mod, dplyr::all_of(feature_col)) %>%
      ggpubr::ggtexttable(rows = c(),
                          cols = c("", tail(colnames(.), n = -1)))
  }
  p <- ggpubr::ggarrange(p_barplot, tbl,
                         labels = c("A", "B"),
                         nrow = 2)
  output <- paste0(opts$options$output, "/", seq_id, ".", opts$options$device)
  ggsave(output, p, width = 20, height = 10)
  saveRDS(p, paste0(opts$options$output, "/", seq_id, ".rds"))
  print(p)
}
