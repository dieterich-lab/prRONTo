library(magrittr)
library(ggplot2)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output"),
  optparse::make_option(c("-s", "--seqid",
                        type = "character",
                        help = "Filter by seqid"))
)

args <- c("--output=mods.pdf",
         "output/data/mods.tsv")
opts <- debug_opts(option_list, args)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t")
df$mod <- as.factor(df$mod)
if (!is.null(opts$options$seqid)) {
  df <- df %>% dplyr::filter(seq_id == opts$options$seqid)
}

p <- df %>% dplyr::group_by(mod) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  ggplot(aes(x = mod, y = count, fill = mod)) +
  geom_col() +
  labs(x = "modification", fill = "") +
  geom_text(aes(label = count)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom")

save_plot(p, opts$options$output, "pdf")
