library(ggplot2)
library(magrittr)

option_list <- list(
  optparse::make_option(c("-f", "--format"),
                        type = "character",
                        default = "pdf_document",
                        help = "Output format"),
  optparse::make_option(c("-p", "--params"),
                        type = "character",
                        help = "YAML file with params"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

params <- yaml::yaml.load_file(opts$options$params)
rmarkdown::render(opts$args,
                  output_format = opts$options$format,
                  params = params)
