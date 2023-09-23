library(magrittr)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

option_list <- list(
  optparse::make_option(c("-a", "--analysis"),
                        type = "character",
                        help = "Analysis type"),
  optparse::make_option(c("-c", "--comparison"),
                        type = "character",
                        help = "Conditions compared"),
  optparse::make_option(c("-f", "--features"),
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
  positional_arguments = TRUE,
  #args=c("--analysis=analysis",
  #       "--comparison=cond1_vs_cond2",
  #       "--features=M,MDI,M_D_I,M5DI,M5D5I5,M_DI,M5_D5_I5,M_M_M5",
  #       "--output=tmp/test.pdf",
  #       "--parameters=preprocessed",
  #       "output/results/analysis/jacusa2/preprocessed/cond1_vs_cond2.out")
)

stopifnot(!is.null(opts$options$output))
stopifnot(opts$options$analysis %in% c("analysis", "downsampling", "mixing"))
stopifnot(!is.null(opts$options$comparison))
stopifnot(nzchar(opts$options$parameters))
stopifnot(nzchar(opts$options$features))
stopifnot(!is.null(opts$args))

add_feature <- function(result, feature, context) {
  feature_col <- paste0("feature_", feature)
  sites <- result[, feature_col]
  regions <- IRanges::resize(GenomicRanges::GRanges(sites), context, fix = "center") %>%
    unique()
  GenomicRanges::mcols(regions) <- NULL

  col <- paste0("feature_", feature, context)
  sites <- plyranges::join_overlap_inner_directed(regions, sites) %>%
    plyranges::group_by(seqnames, start, end, strand) %>%
    plyranges::summarise(
      !!col := sum(!!rlang::ensym(feature_col), na.rm = TRUE)
    ) %>%
    GenomicRanges::GRanges() %>%
    IRanges::resize(1, fix = "center")

  plyranges::join_overlap_inner_directed(result, sites)
}

result <- JACUSA2helper::read_result(opts$args, unpack = TRUE) %>%
  plyranges::select(score, insertion_score, deletion_score, ref)

GenomicRanges::mcols(result, level = "within")[, "feature_M"] <-
  GenomicRanges::mcols(result, level = "within")[, "score"]
GenomicRanges::mcols(result, level = "within")[, "feature_D"] <-
  GenomicRanges::mcols(result, level = "within")[, "deletion_score"]
GenomicRanges::mcols(result, level = "within")[, "feature_I"] <-
  GenomicRanges::mcols(result, level = "within")[, "insertion_score"]

opt_cols <- c()
if (!is.null(opts$options$comparison)) {
  result <- result %>%
    plyranges::mutate(
      comparison = opts$options$comparison)
  opt_cols <- c(opt_cols, "comparison")
}
if (!is.null(opts$options$analysis)) {
  result <- result %>%
    plyranges::mutate(
      analysis = opts$options$analysis)
  opt_cols <- c(opt_cols, "analysis")
}
if (!is.null(opts$options$parameters)) {
  result <- result %>%
    plyranges::mutate(
      parameters = opts$options$parameters)
  opt_cols <- c(opt_cols, "parameters")
}

result <- result %>%
  plyranges::mutate(
    feature_M = replace(feature_M, is.na(feature_M) | feature_M < 0, 0),
    feature_D = replace(feature_D, is.na(feature_D) | feature_D < 0, 0),
    feature_I = replace(feature_I, is.na(feature_I) | feature_I < 0, 0)
  )

features <- strsplit(opts$options$feature, ",") %>%
  unlist(use.names = FALSE)
for (feature in features) {
  feature_col <- paste0("feature_", feature)
  if (feature_col %in% names(GenomicRanges::mcols(result))) {
    next
  }

  sub_features <- strsplit(feature, "_") %>%
    unlist(use.names = FALSE)
  for (sub_feature in sub_features) {
    sub_feature_col <- paste0("feature_", sub_feature)
    if (sub_feature_col %in% names(GenomicRanges::mcols(result))) {
      next
    }

    subsub_features <- stringr::str_extract_all(sub_feature, "([MDI]{1})([0-9]?)") %>%
      unlist(use.names = FALSE)
    for (subsub_feature in subsub_features) {
      subsub_feature_col <- paste0("feature_", subsub_feature)
      if (!subsub_feature_col %in% names(GenomicRanges::mcols(result))) {
        m <- stringr::str_match_all(subsub_feature, "([MDI]{1})([0-9]?)") %>%
          unlist(use.names = FALSE)
        context <- m[3]
        if (context == "") {
          context <- 1
        } else {
          context <- as.numeric(context)
        }
        result <- add_feature(result, m[2], context)
      }
    }
    GenomicRanges::mcols(result)[, sub_feature_col] <- rowSums(
      as.data.frame(GenomicRanges::mcols(result)[, paste0("feature_", subsub_features)]))
  }
}

df <- result %>%
  as.data.frame() %>%
  dplyr::select(
    seqnames,
    end,
    strand,
    ref,
    dplyr::starts_with("feature_"), # FIXME remove not needed
    dplyr::all_of(opt_cols)
  ) %>%
  dplyr::rename(pos = end)

write.table(df, opts$options$output,
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
