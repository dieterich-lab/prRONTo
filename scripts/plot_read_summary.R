library(magrittr)
library(ggplot2)


option_list <- list(
  optparse::make_option(c("-d", "--device"),
                        type = "character",
                        default = "pdf",
                        help = "Plot type"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE#,
  #args = c("--output=tmp",
  #         "output/results/read_summary")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$args))

df <- read.table(opts$args, header = TRUE, sep = "\t", comment.char = "", check.names = FALSE) %>%
  dplyr::mutate(condition = dplyr::case_match(condition, 1 ~ "cond1", 2 ~ "cond2"))
colnames(df) <- gsub("#", "", colnames(df))
df <- df %>%
  dplyr::mutate(analysis = dplyr::case_when(
    parameters %in% c("raw", "preprocessed") ~ "data",
    grepl("^seed~.+_reads~[0-9]+$", parameters) ~ "downsampling",
    .default = "unknown"))

# TODO
# (mixing)

downsampling <- stringr::str_match(df$parameters, "seed~(.+)_reads~([0-9]+)") %>%
  as.data.frame()
colnames(downsampling) <- c("parameters", "seed", "reads")
downsampling <- downsampling %>%
  dplyr::select(-parameters)
df <- dplyr::bind_cols(df, downsampling)
df$parameters <- forcats::fct_relevel(df$parameters, "raw", "preprocessed")

plot_data <- function(df) {
  p <- df %>%
    ggplot(aes(x = interaction(replicate, condition),
               y = numreads,
               fill = parameters)) +
      geom_col(position = position_dodge(width = NULL)) +
      geom_text(aes(label = numreads), position = position_dodge(width = 0.9)) +
      labs(x = "conditions and replicates", y = "reads", fill = "reads", colour = "") +
      scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
      ggtitle("Observed") +
      theme_bw() +
      theme(legend.position = "bottom")

  p
}

plot_seqids <- function(df) {
  p <- df %>%
    ggplot(aes(x = interaction(replicate, condition),
               y = numreads,
               fill = parameters)) +
      geom_col(position = position_dodge(width = NULL)) +
      geom_text(aes(label = numreads), position = position_dodge(width = 0.9)) +
      labs(x = "conditions and replicates", y = "reads", fill = "reads", colour = "") +
      scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
      ggtitle("Sequence IDs") +
      theme_bw() +
      theme(legend.position = "bottom") +
      facet_wrap(rname ~ .)

  p
}

plot_downsampling <- function(df) {
  p <- df %>%
    ggplot(aes(x = interaction(replicate, condition),
               y = numreads)) +
      geom_violin() +
      geom_point(aes(colour = seed)) +
      geom_hline(aes(yintercept = reads, linetype = "dashed")) +
      facet_wrap(~ reads, ncol = 1, labeller = label_both) +
      labs(x = "conditions and replicates", y = "reads", colour = "seed") +
      guides(linetype = "none") +
      scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
      ggtitle("Downsampling") +
      theme_bw() +
      theme(legend.position = "bottom")

  p
}

plot <- function(df, output) {
  p_data <- df %>%
    dplyr::filter(analysis == "data") %>%
    dplyr::group_by(parameters, condition, replicate) %>%
    dplyr::summarise(numreads = sum(numreads)) %>%
    dplyr::ungroup() %>%
    plot_data()

  p_seqids <- df %>%
    dplyr::filter(analysis == "data") %>%
    dplyr::group_by(rname, parameters, condition, replicate) %>%
    dplyr::summarise(numreads = sum(numreads)) %>%
    dplyr::ungroup() %>%
    plot_seqids()

  p_downsampling <- df %>%
    dplyr::filter(analysis == "downsampling") %>%
    dplyr::group_by(seed, reads, parameters, condition, replicate) %>%
    dplyr::summarise(numreads = sum(numreads)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(reads = as.numeric(reads)) %>%
    plot_downsampling()

  p <- ggpubr::ggarrange(p_data, p_seqids, p_downsampling,
                         labels = c("A", "B", "C"),
                         ncol = 1)

  ggsave(output, p, width = 10, height = 20)
}

dir.create(opts$options$output, showWarnings = FALSE)

plot(df, paste0(opts$options$output, "/total.", opts$options$device))
l <- split(df, df$rname)
for (seq_id in names(l)) {
  plot(l[[seq_id]], paste0(opts$options$output, "/", seq_id, ".", opts$options$device))
}
