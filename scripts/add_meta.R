library(magrittr)

# TODO 0, 1 index

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
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$options$output))
stopifnot(opts$options$mods | opts$options$fasta)
stopifnot(!opts$options$fasta & opts$options$context)
stopifnot(!is.null(opts$args))

add_ref_context <- function(result, fasta_fname, context) {
  fasta <- Biostrings::readDNAStringSet(opts$options$fasta_fname)

  GenomeInfoDb::seqlevels(result) <- GenomeInfoDb::seqlevels(fasta)
  GenomicRanges::seqinfo(result) <- GenomicRanges::seqinfo(fasta)

  # TODO test parameters
  result %>%
    IRanges::resize(width = context, fix = "center") %>%
    IRanges::shift(1) %>%
    plyranges::filter(start > 0 & end < GenomeInfoDb::seqlengths(.)[as.character(GenomeInfoDb::seqnames(.))]) %>%
    plyranges::mutate(ref_context = BSgenome::getSeq(fasta, .) %>%
                      as.character()) %>%
    IRanges::shift(-1) %>%
    IRanges::resize(width = IRanges::width(.) + 1)
}

add_mods <- function(result, mods_fname) {
  mods <- data.table::fread(mods_fname, header = TRUE) %>%
    as.data.frame() %>%
    dplyr::mutate(end = start) %>%
    GenomicRanges::GRanges()

  suppressWarnings({
    result %>%
      plyranges::join_overlap_left(mods) })
}

opt_cols <- c()

result <- read.table(opts$args, sep = "\t", header = TRUE) %>%
  GenomicRanges::GenomicRanges(result) %>%
  IRanges::shift(-1)

if (opts$options$mods) {
  opt_cols <- c(opt_cols, "mod")

  result <- add_mods(result, opts$options$mods) %>%
    plyranges::mutate(
      mod = replace(mod, is.na(mod), "")
    )
}

if (opts$options$fasta) {
  opt_cols <- c(opt_cols, "ref_context")

  result <- add_ref_context(result, opts$options$fasta, opts$options$context)
}

df <- result %>%
  as.data.frame()

write.table(df, opts$options$output,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)
