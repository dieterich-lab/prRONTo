library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output"),
  optparse::make_option(c("-x", "--xvalue"),
                        type = "character",
                        help = "Name of column for x value"),
  optparse::make_option(c("-y", "--yvalue"),
                        type = "character",
                        help = "Name of column for y value")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$xvalue))
stopifnot(!is.null(opts$options$yvalue))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t") %>%
  dplyr::mutate(condition = dplyr::case_match(condition, 1 ~ "condition 1", 2 ~ "condition 2"),
                replicate = factor(replicate),
                bam_type = factor(bam_type, levels = c("raw", "preprocessed"), ordered = TRUE))

# TODO finish plot add ecdf
p <- df %>%
  ggplot(aes(x = .data[[opts$options$xvalue]],
             y = .data[[opts$options$yvalue]],
             group = interaction(replicate, condition),
             colour = replicate)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_grid(condition ~ bam_type)

save_plot(p,
          opts$options$output,
          suffix = stringr::str_extract(opts$options$output, "([^.]+)$"))
