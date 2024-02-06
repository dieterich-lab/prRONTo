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
args <- c(
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
  strsplit(":") %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::rename(neighbors = "V1", contamination = "V2") %>%
  transform(neighbors = as.numeric(neighbors), contamination = as.numeric(contamination))

cols <- c("is_outlier", "is_modified", "is_target", "is_in_neighborhood", "mod")
outlier_col <- paste0("lof_outlier_", features)
df$is_outlier <- 0
df[df[, outlier_col] == 1, "is_outlier"] <- 1
df$short_parameters <- df$parameters %>%
  gsub("reads", "r", .) %>%
  gsub("seed", "s", .)
short_parameters <- rev(unique(df$short_parameters))
filtered <- df %>%
  dplyr::filter(lof_neighbors == lof_params$neighbors,
                lof_contamination == lof_params$contamination)

analysis <- data.frame(set = filtered$short_parameters) %>%
  dplyr::mutate(label = dplyr::case_match(set,
                                          "preprocessed" ~ "full data",
                                          .default = "downsampled data"))

filtered <- filtered %>%
  dplyr::filter(is_outlier == 1) %>%
  dplyr::select(seqnames, pos, strand, short_parameters, dplyr::all_of(cols)) %>%
  dplyr::mutate(value = TRUE) %>%
  tidyr::pivot_wider(id_cols = c(seqnames, pos, strand, is_modified, mod, is_target, is_in_neighborhood),
                     names_from = short_parameters,
                     values_from = value,
                     values_fill = FALSE)
  filtered$outlier_type <- ifelse(filtered$is_modified, "known modification",
                                  ifelse(filtered$is_in_neighborhood != "", "modification in vicinity",
                                 "unknown"))
  filtered$label <- ifelse(filtered$is_modified,
                           paste0(filtered$pos, "\n(", filtered$mod, ")"),
                           ifelse(filtered$is_in_neighborhood != "",
                                  paste0(filtered$pos, "\n{", gsub(",", "\n", filtered$is_in_neighborhood), "}"),
                                  filtered$pos))

  p <- upset(filtered, short_parameters,
             name = "Identified in data",
             base_annotations = list(
               "Identified outlier" = intersection_size(counts = FALSE, mapping = aes(fill = outlier_type)) +
               geom_text(
                 mapping = aes(label = label),
                 position = position_stack(vjust = 0.5),
                 na.rm = TRUE) +
               labs(fill = "outlier type") +
               scale_fill_manual(values = c("known modification" = "#66c2a5",
                                            "unknown" = "#fc8d62",
                                            "modification in vicinity" = "#8da0cb")) +
               scale_y_continuous(breaks = integer_breaks()) +
               theme(legend.position = "bottom")
             ),
             width = 0.1,
             set_sizes = FALSE,
             themes = upset_default_themes(legend.position = "bottom"),
             stripes = upset_stripes(mapping = aes(color = label),
                                     geom = geom_segment(size = 4),
                                     colors = c("full data" = "lightpink",
                                                "downsampled data" = "snow3"),
                                     data = analysis)
) +
  ylab("Data (BAM files)") +
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 1))

save_plot(p, opts$options$output, "pdf")
