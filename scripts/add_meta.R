library(magrittr)
source(paste0(Sys.getenv("PRONTO_DIR"), "/scripts/pronto_lib.R"))

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
  optparse::make_option(c("-n", "--neighborhood"),
                        type = "integer",
                        help = "Up/downstream neighborhood to a known modification"),
  optparse::make_option(c("-t", "--targets"),
                        type = "character",
                        help = "','-separated list of targets"),
  optparse::make_option(c("-s", "--fasta"),
                        type = "character",
                        help = "Reference sequence FASTA")
)
args <- c(
  "--output=test/results/analysis/jacusa2/preprocessed/cond1_vs_cond2_meta.tsv",
  "--mods=test/data/mods.tsv",
  "--fasta=test/data/ref.fasta",
  "--context=5",
  "--neighborhood=2",
  "test/results/downsampling/jacusa2/seed~3_reads~200/cond1_vs_cond2_scores.tsv")
opts <- debug_opts(option_list, args)

stopifnot(!is.null(opts$options$output))
stopifnot(!is.null(opts$options$mods) || !is.null(opts$options$fasta))
stopifnot(!(is.null(opts$options$fasta) && !is.null(opts$options$context)))
stopifnot(!(is.null(opts$options$mods) && !is.null(opts$options$neighborhood)))
stopifnot(!is.null(opts$args))

add_targets <- function(result, targets) {
  targets <- strsplit(targets, ",") %>%
    unlist() %>%
    GenomicRanges::GRanges()


  suppressWarnings({
    result %>%
      plyranges::join_overlap(targets)})
}

# TODO test
add_neighborhood <- function(result, mods, neighborhood) {
  region <- result
  GenomicRanges::mcols(region) <- NULL

  region <- region %>%
    plyranges::mutate(site = as.character(result), site_pos = start) %>%
    IRanges::resize(width = 2 * neighborhood + 1, fix = "center") %>%
    plyranges::find_overlaps(mods) %>%
    plyranges::filter(site_pos != mod_pos) %>%
    plyranges::mutate(mod = paste0(mod_pos, ":", mod)) %>%
    plyranges::select(site, mod) %>%
    GenomicRanges::mcols() %>%
    as.data.frame() %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(mods = unique(mod) %>% paste0(collapse = ",")) %>%
    dplyr::select(site, mods) %>%
    dplyr::distinct() %>%
    dplyr::rename(is_in_neighborhood = mods) %>%
    tidyr::separate(col = site, into = c("seqnames", "start"), sep = ":") %>%
    dplyr::mutate(end = start) %>%
    GenomicRanges::GRanges()

  plyranges::join_overlap_left(result, region)
}

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

get_mods <- function(mods_fname) {
  data.table::fread(mods_fname, header = TRUE) %>%
    as.data.frame() %>%
    dplyr::rename(seqnames = seq_id, start = pos) %>%
    dplyr::mutate(start = start + 1, end = start, mod_pos = start) %>%
    GenomicRanges::GRanges()
}

add_mods <- function(result, mods) {
  suppressWarnings({
    result %>%
      plyranges::join_overlap_left(mods)})
}

opt_cols <- c()

result <- read.table(opts$args, sep = "\t", header = TRUE) %>%
  dplyr::mutate(start = pos, end = pos) %>%
  GenomicRanges::GRanges()

GenomicRanges::mcols(result)["is_modified"] <- 0
if (!is.null(opts$options$mods)) {
  opt_cols <- c(opt_cols, "mod")

  mods <- get_mods(opts$options$mods)
  result <- add_mods(result, mods) %>%
    plyranges::mutate(
      mod = replace(mod, is.na(mod), "*"),
      is_modified = dplyr::case_match(mod, "*" ~ 0, .default = 1)
    )
  if (!is.null(opts$options$neighborhood)) {
    opt_cols <- c(opt_cols, "is_in_neighborhood")
    result <- add_neighborhood(result, mods, opts$options$neighborhood)
  } else {
    GenomicRanges::mcols(result)["is_in_neighborhood"] <- NA
  }
} else {
  GenomicRanges::mcols(result)["is_in_neighborhood"] <- NA
}

if (!is.null(opts$options$fasta)) {
  opt_cols <- c(opt_cols, "ref_context")

  result <- add_ref_context(result, opts$options$fasta, opts$options$context)
}

GenomicRanges::mcols(result)["is_target"] <- 0
if (!is.null(opts$options$targets)) {
  results <- add_targets(targets)
}

df <- result %>%
  as.data.frame() %>%
  dplyr::select(-c(start, end, width)) %>%
  dplyr::relocate(ref_context, .after = ref)

write.table(df, opts$options$output,
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
