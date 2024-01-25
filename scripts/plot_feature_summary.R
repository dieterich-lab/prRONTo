library(ggplot2)
library(ComplexUpset)
library(magrittr)

# FIXME remove
#Sys.setenv(PRONTO_DIR = "~/mnt/beegfs/homes/mpiechotta/git/prRONTo/",
#           PRONTO_DEBUG = "TRUE")
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-l", "--lof_params"),
                        type = "character",
                        help = "','-separted list of neighbors:contamination"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
args = c("--lof_params=20:0.001",
         "--output=tmp/tmp.pdf",
         "~/mnt/beegfs/homes/mpiechotta/git/prRONTo/example/human/output_test/results/merged_lof.tsv")
opts = debug_opts(option_list, args)
stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))
# TODO add sensible stopifnot

df <- read_merged_lof(opts$args) %>%
  dplyr::mutate(lof_params = paste0(lof_neighbors, ":", lof_contamination))
df$lof <- paste0("n=", df$lof_neighbors, ", c=", df$lof_contamination)
df$short_parameters <- df$parameters %>%
  gsub("reads", "r", .) %>%
  gsub("seed", "s", .)

# filter by lof_parameters
lof_opts <- strsplit(opts$options$lof_params, ",") %>%
  unlist()

features <- guess_features(df)
outlier_col <- paste0("lof_outlier_", features)
cols <- c(outlier_col, "is_modified", "is_target", "is_in_neighborhood", "mod")
parameters <- rev(unique(df$parameters))
filtered <- df %>%
  dplyr::filter(lof_params %in% lof_opts) %>%
  dplyr::select(seqnames, pos, strand, parameters, dplyr::all_of(cols)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(outlier_col), names_to = "feature", values_to = "is_outlier") %>%
  dplyr::mutate(feature = gsub("outlier_col_", "", feature)) %>%
  dplyr::filter(is_outlier == 1) %>%
  tidyr::pivot_wider(id_cols = c(seqnames, pos, strand, is_modified, mod, is_target, is_in_neighborhood),
                     names_from = parameters,
                     values_from = is_outlier,
                     values_fn = function(x) { as.logical(sum(x)) },
                     values_fill = FALSE)
filtered$outlier_type <- ifelse(filtered$is_modified, "known modification",
                                ifelse(filtered$is_in_neighborhood != "", "modification in vicinity",
                               "unknown"))
filtered$label <- ifelse(filtered$is_modified,
                         paste0(filtered$pos, "\n(", filtered$mod, ")"),
                         ifelse(filtered$is_in_neighborhood != "",
                                paste0(filtered$pos, "\n{", gsub(",", "\n", filtered$is_in_neighborhood), "}"),
                                filtered$pos))

p <- upset(filtered, parameters,
           name = "runs",
           base_annotations = list(
             "Intersection size" = intersection_size(counts = FALSE, mapping = aes(fill = outlier_type)) +
             geom_text(
               mapping = aes(label = label),
               position = position_stack(vjust = 0.5),
               na.rm = TRUE) +
             labs(fill = "outlier type") +
             scale_fill_manual(values = c("known modification" = "#66c2a5",
                                          "unknown" = "#fc8d62",
                                          "modification in vicinity" = "#8da0cb")) +
             scale_y_continuous(breaks = integer_breaks())
           ),
           width = 0.1,
           set_sizes = FALSE) +
      theme(legend.position = "bottom")
save_plot(p, opts$options$output, "pdf")
