library(magrittr)
library(ggplot2)

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
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE,
  args = c("--output=tmp",
           "output/results/merged_lof.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t")

downsampling <- stringr::str_match(df$parameters, "seed~(.+)_reads~([0-9]+)") %>%
  as.data.frame()
colnames(downsampling) <- c("parameters", "seed", "reads")
downsampling <- downsampling %>%
  dplyr::select(-parameters)
df <- dplyr::bind_cols(df, downsampling)

p <- df %>%
  dplyr::filter(analysis == "analysis") %>%
  dplyr::select(seqnames, pos, strand, dplyr::starts_with("feature_")) %>%
  tidyr::pivot_longer(c("feature_M", "feature_MDI"), names_to = "feature", values_to = "score") %>%
  dplyr::mutate(feature = gsub("feature_", "", feature)) %>%
  ggplot(aes(x = score, colour = feature)) +
    geom_density() +
    labs(x = "score") +
    theme_bw() +
    theme(legend.position = "bottom")

# TODO
# downsampling
# seed, reads
# density and confidence intervals

ggsave(opts$options$output, p)
