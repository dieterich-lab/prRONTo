library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-1", "--cond1"),
                        type = "character",
                        default = "condition 1",
                        help = "label for condition"),
  optparse::make_option(c("-2", "--cond2"),
                        type = "character",
                        default = "condition 2",
                        help = "label for condition"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output"),
  optparse::make_option(c("-l", "--label"),
                        type = "character",
                        help = "Label for column"),
  optparse::make_option(c("-n", "--normalize"),
                        type = "character",
                        help = "Column to normalize"),
  optparse::make_option(c("-c", "--column"),
                        type = "character",
                        help = "Columni to use")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$label))
stopifnot(!is.null(opts$options$column))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t", check.names = FALSE) %>%
  dplyr::mutate(condition = dplyr::case_match(condition, 1 ~ opts$options$cond1, 2 ~ opts$options$cond2),
                replicate = factor(replicate),
                bam_type = factor(bam_type, levels = c("raw", "preprocessed"), ordered = TRUE))

col <- opts$options$column
if (col == "reads mapped") {
  df <- df %>%
    dplyr::filter(bam_type != "preprocessed") %>%
    dplyr::mutate(bam_type = factor(bam_type, levels = c("raw"), ordered = TRUE))
}

if (!is.null(opts$options$normalize)) {
  new_col <- paste0("norm_", opts$options$column)
  df[, new_col] <- round(df[, opts$options$column] / df[, opts$options$normalize], digits = 4)
  col <- new_col
}

p <- df %>%
  ggplot(aes(x = replicate,
             y = .data[[col]],
             group = interaction(replicate, condition),
             fill = replicate,
             colour = replicate)) +
    geom_col() +
    labs(y = opts$options$label) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_grid(bam_type ~ condition)

if (!is.null(opts$options$normalize)) {
  p <- p + scale_y_continuous(labels = scales::percent)
}


save_plot(p,
          opts$options$output,
          suffix = stringr::str_extract(opts$options$output, "([^.]+)$"))
