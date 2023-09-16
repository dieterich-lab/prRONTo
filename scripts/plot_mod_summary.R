library(magrittr)
library(ggplot2)

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)

opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE#,
  #args = c("--output=mods.pdf",
  #         "output/data/mods.tsv")
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
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

# TODO
# covered in bams

ggsave(opts$options$output, p)
saveRDS(p, paste0(gsub(".pdf$", "", opts$options$output), ".rds"))
