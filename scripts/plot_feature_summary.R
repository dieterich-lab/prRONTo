library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE#,
  #args = c("--output=tmp/test.pdf",
  #         "output/results/merged_lof.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read_merged_lof(opts$args)

filtered <- df %>%
  #dplyr::filter(analysis == "analysis") %>%
  dplyr::select(seqnames, pos, strand, dplyr::starts_with("lof_outlier_"), mod, seed, reads) %>%
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
  tidyr::unite(col = "lof_outlier", dplyr::all_of(lof_outliers), sep = ",", na.rm = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(lof_outlier = strsplit(lof_outlier, ","))

p <- filtered %>%
  ggplot(aes(x = lof_outlier, fill = is_mod)) +
  geom_bar() +
  labs(x = "feature score", fill = "outlier") +
  ggupset::scale_x_upset() +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot(p, opts$options$output)
