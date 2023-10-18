library(ggplot2)
library(ComplexUpset)
library(magrittr)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-f", "--features"),
                        type = "character",
                        help = "','-separted list of features, otherwise guessed from column names"),
  optparse::make_option(c("-a", "--analysis"),
                        type = "character",
                        help = "','-separted list of original and,or downsampling"),
  optparse::make_option(c("-l", "--lof_params"),
                        type = "character",
                        help = "','-separted list of neighbors:contamination"),
  optparse::make_option(c("-r", "--regions"),
                        type = "character",
                        help = "','-separted list of regions"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
args = c(
         "--features=M",
         "--lof=20:0.001",
         "--output=tmp",
         "output/results/merged_lof.tsv")
opts <- debug_opts(option_list, args)
stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))
# TODO add sensible stopifnot

df <- read_merged_lof(opts$args)

regions <- NULL
if (!is.null(opts$options$regions)) {
  regions <- unlist(strsplit(opts$options$regions, ","))
  df <- df %>%
    dplyr::filter(seqnames %in% regions)
}
features <- unlist(strsplit(opts$options$features, ","))

lof_params <- strsplit(opts$options$lof_params, ",") %>%
  unlist() %>%
  strsplit(":") %>% do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::rename(neighbors = "V1", contamination = "V2") %>%
  transform(neighbors = as.numeric(neighbors), contamination = as.numeric(contamination))

if (length(features) == 1) {
  cols <- c("is_outlier", "is_modified", "is_target", "is_in_neighborhood")
  outlier_col <- paste0("lof_outlier_", features)
  df$is_outlier <- 0
  df[df[, outlier_col] == 1, "is_outlier"] <- 1
  parameters <- rev(unique(df$parameters))
  filtered <- df %>%
    dplyr::filter(is_outlier == 1, lof_neighbors == lof_params$neighbors, lof_contamination == lof_params$contamination) %>%
    dplyr::select(seqnames, pos, strand, parameters, dplyr::all_of(cols)) %>%
    dplyr::mutate(value = TRUE) %>%
    tidyr::pivot_wider(id_cols = c(seqnames, pos, strand, is_modified, is_target, is_in_neighborhood),
                       names_from = parameters,
                       values_from = value,
                       values_fill = FALSE)
    filtered$outlier_type <- ifelse(filtered$is_modified, "known modification",
                            ifelse(filtered$is_in_neighborhood, "modification in vicinity",
                                   "unknown"))

    p <- upset(filtered, parameters,
               name = "runs",
               base_annotations = list(
                 "Intersection size" = intersection_size(counts = FALSE, mapping = aes(fill = outlier_type)) +
                   labs(fill = "outlier type") +
                   scale_fill_manual(values = c("known modification" = "#66c2a5",
                                                "unknown" = "#fc8d62",
                                                "modification in vicinity" = "#8da0cb"))
               ),
               width = 0.1,
               set_sizes = FALSE
    )
}
# TODO -> feature, LOFs

save_plot(p, opts$options$output, "pdf")
