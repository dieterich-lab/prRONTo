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
  args = c("--output=mods.pdf",
           "output/data/mods.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t")

p <- df %>% dplyr::group_by(mod) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  ggplot(aes(x = mod, y = count, fill = mod)) +
  geom_col() +
  labs(x = "modification", fill = "") +
  geom_text(aes(label = count)) +
  theme_bw() +
  theme(legend.position = "bottom")

# TODO
# covered in bams

ggsave(opts$options$output, p)
