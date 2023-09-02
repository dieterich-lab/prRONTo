library(magrittr)
options(error = traceback)

option_list <- list(
  optparse::make_option(c("-a", "--analysis"),
                        type = "character",
                        help = "Analysis type"),
  optparse::make_option(c("-c", "--comparison"),
                        type = "character",
                        help = "Conditions compared"),
  optparse::make_option(c("-s", "--scores"),
                        type = "character",
                        help = "Features to calculate"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output"),
  optparse::make_option(c("-p", "--parameters"),
                        type = "character",
                        help = "Parameters")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(opts$options$analysis %in% c("analysis", "downsampling", "mixing"))
stopifnot(opts$options$comparison %in% c("cond1_vs_cond2", "cond1_vs_cond3", "cond2_vs_cond3"))
stopifnot(nzchar(opts$options$parameters))
stopifnot(nzchar(opts$options$scores))
stopifnot(!is.null(opts$args))

## move to JACUSA2helper

add_score <- function(result, score, context, col) {
  sites <- result[, score]
  regions <- IRanges::resize(GenomicRanges::GRanges(sites), context, fix = "center") %>%
    unique()
  GenomicRanges::mcols(regions) <- NULL

  sites <- plyranges::join_overlap_inner_directed(regions, sites) %>%
    plyranges::group_by(seqnames, start, end, strand) %>%
    plyranges::summarise(
      !!col := sum(!!rlang::ensym(score), na.rm = TRUE)
    ) %>%
    GenomicRanges::GRanges() %>%
    IRanges::resize(1, fix = "center")

  plyranges::join_overlap_inner_directed(result, sites)
}

add_scores <- function(result, scores) {
  for (new_score in names(scores)) {
    if (!new_score %in% names(new_score)) {
      score_context <- scores[[new_score]]
      for (j in 1:nrow(score_context)) {
        score <- score_context$score[j]
        context <- score_context$context[j]
        col <- score_context$col[j]
        if (context > 1) {
          result <- add_score(result, score, context, col)
        }
      }
      new_values <- rowSums(as.data.frame(GenomicRanges::mcols(result[, score_context$col])))
      result <- result %>%
        plyranges::mutate(!!new_score := new_values)
    }
  }

  result
}

result <- JACUSA2helper::read_result(opts$args, unpack = TRUE) %>%
  IRanges::shift(-1) %>% # TODO JACUSA2 stores 0-indexed
  plyranges::select(score, insertion_score, deletion_score, ref)

GenomicRanges::mcols(result, level = "within")[, "M"] <-
  GenomicRanges::mcols(result, level = "within")[, "score"]
GenomicRanges::mcols(result, level = "within")[, "D"] <-
  GenomicRanges::mcols(result, level = "within")[, "deletion_score"]
GenomicRanges::mcols(result, level = "within")[, "I"] <-
  GenomicRanges::mcols(result, level = "within")[, "insertion_score"]

opt_cols <- c()
if (!is.null(opts$options$comparison)) {
  result <- result %>%
    plyranges::mutate(
      comparison = opts$options$comparison)
  opt_cols <- "comparison"
}
if (!is.null(opts$options$analysis)) {
  result <- result %>%
    plyranges::mutate(
      analysis = opts$options$analysis)
  opt_cols <- "analysis"
}
if (!is.null(opts$options$parameters)) {
  result <- result %>%
    plyranges::mutate(
      parameters = opts$options$parameters)
  opt_cols <- "parameters"
}

result <- result %>%
  plyranges::mutate(
    M = replace(M, is.na(M) | M < 0, 0),
    D = replace(D, is.na(D) | D < 0, 0),
    I = replace(I, is.na(I) | I < 0, 0)
  )

# parse scores str to list of score and context data frames
score_cols <- strsplit(opts$options$scores, ",") %>%
  unlist()

scores <- score_cols %>%
  stringr::str_match_all("([MDI]{1})([0-9]?)") %>%
    lapply(function(m) {
             df <- data.frame(score = m[, 2], context = as.numeric(m[, 3]))
             df[is.na(df$context), "context"] <- 1
             df$col <- paste0(df$score, df$context)
             df$col <- gsub("1", "", df$col)

             return(df)
  })
names(scores) <- score_cols

result <- add_scores(result, scores)

df <- result %>%
  as.data.frame() %>%
  dplyr::select(
    seqnames,
    end,
    strand,
    dplyr::all_of(score_cols),
    dplyr::all_of(opt_cols),
    -c(score, deletion_score, insertion_score)
  ) %>%
  dplyr::rename(pos = end)

write.table(df, opts$options$output,
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
