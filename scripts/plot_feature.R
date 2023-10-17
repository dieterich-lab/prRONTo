# TODO
# * outlier
# * vicintiy
# * targets, missed modifications
# * FIXME no outliers
# * knownmodifcation

library(ggplot2)
library(magrittr)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))


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
         "output/results/analysis/jacusa2/preprocessed/lof/neighbors~20_contamination~0.002/cond1_vs_cond2.tsv")
opts <- debug_opts(option_list, args)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$feature))
stopifnot(!is.null(opts$args))

feature <- opts$options$feature
feature_col <- paste0("feature_", feature)
lof_outlier_col <- paste0("lof_outlier_", feature)

result <- read.table(opts$args, sep = "\t", header = TRUE) %>%
  dplyr::mutate(is_outlier = dplyr::case_when(.data[[lof_outlier_col]] == 1 ~ TRUE, .default = FALSE),
                label = dplyr::case_when(is_modified == 1 & is_outlier == 1 ~ paste(pos,"(", mod, ")"),
                                         is_modified == 0 & is_outlier == 1 ~ as.character(pos)),
                outlier_type = dplyr::case_when(is_modified == 1 & is_outlier == 1 ~ "known modification",
                                                is_modified == 0 & is_outlier == 1 ~ "unknown",
                                                is_in_neighborhood == 1 & is_outlier == 1 ~ "modification in vicinity",
                                                .default = "inlier")) %>%
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
    ggplot(aes(x = pos, y = .data[[feature_col]], colour = outlier_type)) +
      scale_colour_manual(values = c("known modification" = "#66c2a5",
                            "unknown" = "#fc8d62",
                            "modification in vicinity" = "#8da0cb",
                            "inlier" = "#c0c0c0")) +
      geom_hline(aes(yintercept = 0)) +
      geom_segment(aes(x = pos,
                       y = .data[[feature_col]],
                       xend = pos,
                       yend = .data[[feature_col]] - .data[[feature_col]])) +
      geom_point(data = outlier, size = 3) +
      geom_text(data = outlier,
                mapping = aes(label = label),
                hjust = 0, vjust = 0) +
      guides(colour = guide_legend(title = "outlier type")) +
      labs(x = "Position [nt]", y = feature) +
      theme_bw() +
      theme(legend.position = "bottom")

  tbl = "no data"
  # FIXME remove
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
  saveRDS(p_barplot, paste0(opts$options$output, "/", seq_id, "_barplot.rds"))
  write.table(outlier,
              paste0(opts$options$output, "/", seq_id, "_outlier.tsv"),
              quote = FALSE,
              row.names = FALSE,
              sep = "\t")
}
