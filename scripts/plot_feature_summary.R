# FIXME

library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
args = c("--output=tmp/test.pdf",
         "output/results/merged_lof.tsv")
opts <- debug_opts(option_list, args)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read_merged_lof(opts$args)
df$lof <- paste0("n=", df$lof_neighbors, ", c=", df$lof_contamination)
df$short_parameters <- df$parameters %>%
  gsub("reads", "r", .) %>%
  gsub("seed", "s", .)


filtered <- df %>%
  #dplyr::filter(analysis == "analysis") %>%
  dplyr::select(seqnames, pos, strand, dplyr::starts_with("lof_outlier_"), mod, short_parameters, lof) %>%
  dplyr::filter(dplyr::if_any(starts_with("lof_outlier_"), ~ .x == 1)) %>%
  dplyr::mutate(outlier_type = dplyr::case_when(mod != "*" ~ "known modification", .default = "unknown"))


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
  ggplot(aes(x = lof_outlier, fill = outlier_type)) +
  geom_bar() +
  labs(x = "feature", fill = "outlier") +
  scale_fill_manual(values = c("known modification" = "#66c2a5",
                        "unknown" = "#fc8d62",
                        "modification in vicinity" = "#8da0cb")) +
  ggupset::scale_x_upset() +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_grid(short_parameters ~ lof)

save_plot(p, opts$options$output, "pdf")
