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

filtered <- df %>%
  dplyr::filter(analysis == "analysis") %>%
  dplyr::select(seqnames, pos, strand, dplyr::starts_with("lof_outlier_"), mod) %>%
  dplyr::filter(dplyr::if_any(starts_with("lof_outlier_"), ~ .x == 1)) %>%
  dplyr::mutate(is_mod = dplyr::case_when(mod != "*" ~ "known mod", .default = "unknown"))

lof_outliers = grep("^lof_outlier_", colnames(filtered), value = TRUE)
for (lof_outlier in lof_outliers) {
  feature <- gsub("^lof_outlier_", "", lof_outlier)
  filtered[[lof_outlier]] <- dplyr::case_match(filtered[[lof_outlier]],
                                               1 ~ feature,
                                               .default = NULL)
}
filtered <- filtered %>%
  tidyr::unite(col = "lof_outlier", lof_outliers, sep = ",", na.rm = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(lof_outlier = strsplit(lof_outlier, ","))

p <- filtered %>%
  ggplot(aes(x = lof_outlier, fill = is_mod)) +
  geom_bar() +
  labs(x = "feature score", fill = "outlier") +
  scale_x_upset() +
  theme_bw() +
  theme(legend.position = "bottom")

#ggsave(opts$options$output, p)
