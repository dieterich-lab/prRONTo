library(magrittr)

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


CONTEXT <- 5

option_list <- list(
  optparse::make_option(c("-m", "--mods"),
                        type = "character",
                        help = "Modifications"),
  optparse::make_option(c("-o", "--output"),
                        type = "character",
                        help = "Output"),
  optparse::make_option(c("-c", "--context"),
                        type = "numeric",
                        default = CONTEXT,
                        help = paste0("Ref context. Default: ", CONTEXT)),
  optparse::make_option(c("-s", "--fasta"),
                        type = "character",
                        help = "Reference sequence FASTA")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE#,
  #args = c(
  #  "--output=output/results/analysis/jacusa2/preprocessed/cond1_vs_cond2_meta.tsv",
  #  "--mods=output/data/mods.tsv",
  #  "--fasta=output/data/ref.fasta",
  #  "--context=5",
    #"output/results/analysis/jacusa2/preprocessed/cond1_vs_cond2_scores.tsv"),
  #  "output/results/downsampling/jacusa2/seed~3_reads~200/cond1_vs_cond2_scores.tsv")
)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$mods) || !is.null(opts$options$fasta))
# FIXME stopifnot(is.null(opts$options$fasta) && !is.null(opts$options$context))
stopifnot(!is.null(opts$args))

add_ref_context <- function(result, fasta_fname, context) {
  fasta <- Biostrings::readDNAStringSet(fasta_fname)

  GenomeInfoDb::seqlevels(result) <- GenomeInfoDb::seqlevels(fasta)
  GenomicRanges::seqinfo(result) <- GenomicRanges::seqinfo(fasta)

  result %>%
    IRanges::resize(width = context, fix = "center") %>%
    plyranges::filter(start > 0 & end <= GenomeInfoDb::seqlengths(.)[as.character(GenomeInfoDb::seqnames(.))]) %>%
    plyranges::mutate(ref_context = BSgenome::getSeq(fasta, .) %>%
                      as.character()) %>%
    IRanges::resize(width = 1, fix = "center")
}

add_mods <- function(result, mods_fname) {
  mods <- data.table::fread(mods_fname, header = TRUE) %>%
    as.data.frame() %>%
    dplyr::rename(seqnames = seq_id, start = pos) %>%
    dplyr::mutate(end = start) %>%
    GenomicRanges::GRanges() %>%
    IRanges::shift(1)

  suppressWarnings({
    result %>%
      plyranges::join_overlap_left(mods)})
}

opt_cols <- c()

result <- read.table(opts$args, sep = "\t", header = TRUE) %>%
  dplyr::mutate(start = pos, end = pos) %>%
  GenomicRanges::GRanges()

if (!is.null(opts$options$mods)) {
  opt_cols <- c(opt_cols, "mod")

  result <- add_mods(result, opts$options$mods) %>%
    plyranges::mutate(
      mod = replace(mod, is.na(mod), "*")
    )
}

if (!is.null(opts$options$fasta)) {
  opt_cols <- c(opt_cols, "ref_context")

  result <- add_ref_context(result, opts$options$fasta, opts$options$context)
  head(result$ref_context)
}

df <- result %>%
  as.data.frame() %>%
  dplyr::select(-c(start, end, width)) %>%
  dplyr::relocate(ref_context, .after = ref)

write.table(df, opts$options$output,
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
