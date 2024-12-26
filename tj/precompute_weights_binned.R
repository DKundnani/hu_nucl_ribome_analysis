# Get the number of proper nucleotides (A, C, G, T) in each window of the reference genome

weights_binned.save <- function(
  window_width,
  overwrite = FALSE
) {
  purrr::walk(
    const.CHROMS,
    function(chrom) {
      file_out <- fn.weights_binned(
        chrom = chrom,
        window_width = window_width
      )
      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }
      GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(
          start = 1,
          end = utils.get_chrom_sizes(chrom)
        ),
        window_width = window_width
      ) |>
      plyranges::slide_ranges(
        width = window_width,
        step = window_width
      ) |>
      GenomicRanges::granges() |>
      x => plyranges::mutate(
          x,
          count = (
            BSgenome::BSgenomeViews(
              BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
              x
            ) |>
            Biostrings::oligonucleotideFrequency(width = 1) |>
            rowSums()
          )
      ) |>
      utils.write_rds(file_out)
    }
  )
}

weights_binned.do_main <- function(
  window_width_list = c(1e3, 1e4, 1e5, 1e6),
  overwrite = FALSE
) {
  list(window_width = window_width_list) |>
  purrr::pwalk(weights_binned.save, overwrite = overwrite)
}
