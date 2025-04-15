# This file saves the regions on the hg38 genome that are
# either gaps (filled with Ns) or repeat regions determined with repeatMasker.
# This is used to exclude such regions from
# certain analyses. We obtain the gaps by directly searching for the letter N
# in each chromosome. We obtain the repeats by using the masks in
# the BSgenome.Hsapiens.UCSC.hg38.masked package. The resulting files
# should be equivalent to what can be downloaded from the UCSC
# genome browser (https://genome.ucsc.edu/cgi-bin/hgTables).

hg38_annotations.gaps <- function(overwrite = FALSE) {
  file_out <- fn.HG38_GAPS
  if (!overwrite && file.exists(file_out)) {
    utils.log("Already exists:", file_out)
    return(invisible())
  }
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  purrr::map(
    const.CHROMS,
    function(chrom) {
      utils.log("Gaps:", chrom)

      genome[[chrom]] |>
      Biostrings::maskMotif("N") |>
      x => Biostrings::masks(x)[[1]] |>
      x => GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(
          start = IRanges::start(x),
          end = IRanges::end(x)
        )
      ) |>
      GenomicRanges::reduce()
    }
  ) |>
  purrr::reduce(plyranges::bind_ranges) |>
  utils.write_rds(file_out)
}

hg38_annotations.repeats <- function(overwrite = FALSE) {
  file_out <- fn.HG38_REPEATS
  if (!overwrite && file.exists(file_out)) {
    utils.log("Already exists:", file_out)
    return(invisible())
  }
  genome <- BSgenome.Hsapiens.UCSC.hg38.masked::BSgenome.Hsapiens.UCSC.hg38.masked
  purrr::map(
    const.CHROMS,
    function(chrom) {
      utils.log("Repeats:", chrom)

      genome[[chrom]] |>
      Biostrings::masks() |>
      x => x[IRanges::desc(x) == "RepeatMasker"] |>
      x => x[[1]] |>
      x => GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(
          start = IRanges::start(x),
          end = IRanges::end(x)
        )
      ) |>
      GenomicRanges::reduce()
    }
  ) |>
  purrr::reduce(plyranges::bind_ranges) |>
  utils.write_rds(file_out)
}

hg38_annotations.do_main <- function(overwrite = FALSE) {
  hg38_annotations.gaps(overwrite)
  hg38_annotations.repeats(overwrite)
}
