# TODO FINISH
# downsampling
# seed, reads
# density and confidence intervals
# what features how to combine

library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-d", "--device"),
                        type = "character",
                        default = "pdf",
                        help = "Plot type"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE#,
  #args = c("--output=tmp",
  #         "output/results/merged_lof.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read_merged_lof(opts$args)

p <- df %>%
  dplyr::filter(analysis == "analysis") %>%
  dplyr::select(seqnames, pos, strand, dplyr::starts_with("feature_")) %>%
  tidyr::pivot_longer(dplyr::starts_with("feature_"), names_to = "feature", values_to = "score") %>%
  dplyr::mutate(feature = gsub("feature_", "", feature)) %>%
  ggplot(aes(x = score, colour = feature)) +
    geom_density() +
    labs(x = "score") +
    theme_bw() +
    theme(legend.position = "bottom")

save_plot(p, opts$options$output, opts$options$device)
