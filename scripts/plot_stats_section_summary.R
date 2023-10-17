library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output"),
  optparse::make_option(c("-1", "--cond1"),
                        type = "character",
                        default = "condition 1",
                        help = "label for condition"),
  optparse::make_option(c("-2", "--cond2"),
                        type = "character",
                        default = "condition 2",
                        help = "label for condition"),
  optparse::make_option(c("-x", "--xvalue"),
                        type = "character",
                        help = "Name of column for x value"),
  optparse::make_option(c("-X", "--xlabel"),
                        type = "character",
                        help = "Label for x value"),
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
  dplyr::mutate(condition = dplyr::case_match(condition,
                                              1 ~ opts$options$cond1,
                                              2 ~ opts$options$cond2),
                replicate = factor(replicate),
                bam_type = factor(bam_type, levels = c("raw", "preprocessed"), ordered = TRUE))

type = "violin"
p <- NULL
if (type == "line") {
  p <- df %>%
    ggplot(aes(x = .data[[opts$options$xvalue]],
               y = .data[[opts$options$yvalue]],
               colour = replicate)) +
      geom_line() +
      theme_bw() +
      theme(legend.position = "bottom") +
      facet_grid(bam_type ~ condition)

  if (!is.null(opts$options$xlabel)) {
    p <- p +
      labs(x = opts$options$xlabel)
  }
} else if (type == "violin") {
    p <- df %>%
      ggplot(aes(x = replicate,
                 y = .data[[opts$options$xvalue]],
                 weight = .data[[opts$options$yvalue]],
                 colour = replicate)) +
        geom_violin() +
        theme_bw() +
        theme(legend.position = "bottom") +
        facet_grid(bam_type ~ condition)
} else if (type == "cumulative") {
  # TODO
}

save_plot(p,
          opts$options$output,
          suffix = stringr::str_extract(opts$options$output, "([^.]+)$"))
