# Get the number of kmers in each window of the reference genome

kmer_binned.save <- function(
  window_width,
  kmer_size,
  overwrite = FALSE
) {
  purrr::walk(
    const.CHROMS,
    function(chrom) {
      file_out <- fn.kmer_binned(
        chrom = chrom,
        window_width = window_width,
        kmer_size = kmer_size
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
        strand = "+",
        window_width = window_width
      ) |>
      plyranges::slide_ranges(
        width = window_width,
        step = window_width
      ) |>
      GenomicRanges::granges() |>
      x => {
        count_plus <- (
          BSgenome::BSgenomeViews(
            BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
            x
          ) |>
          Biostrings::oligonucleotideFrequency(width = kmer_size)
        )
        count_minus <- count_plus
        colnames(count_minus) <- kmer.reverse_complement(colnames(count_minus))
        count_minus <- count_minus[, colnames(count_plus)]
        
        x_plus <- x
        x_plus$count <- count_plus
        GenomicRanges::strand(x_plus) <- "+"
        x_minus <- x
        x_minus$count <- count_minus
        GenomicRanges::strand(x_minus) <- "-"
        c(x_plus, x_minus)
      } |>
      utils.write_rds(file_out)
    }
  )
}

kmer_binned.do_main <- function(
  window_width_list = c(1e3, 1e4, 1e5, 1e6),
  kmer_size_list = c(1, 2),
  overwrite = FALSE
) {
  tidyr::expand_grid(
    window_width = window_width_list,
    kmer_size = kmer_size_list
  ) |>
  purrr::pwalk(kmer_binned.save, overwrite = overwrite)
}
