# Get ribo data aggregated into bins across the genome

dnaseq_binned.save <- function(window_width, overwrite = FALSE) {
  purrr::pwalk(
    tidyr::expand_grid(
      sample = const.SAMPLES_DNASEQ,
      chrom = const.CHROMS
    ),
    function(sample, chrom) {
      file_out <- fn.dnaseq_binned(
        sample = sample,
        chrom = chrom,
        window_width = window_width,
        assembly = "hg38"
      )
      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }
      granges_dnaseq <- utils.load_dnaseq_rds(sample, chrom)
      if (is.null(granges_dnaseq) || (NROW(granges_dnaseq) == 0)) {
        utils.log("No DNA-seq for", sample, chrom)
        utils.write_rds(NULL, file_out)
        return(invisible())
      }
      granges_binned <- utils.get_ranges_binned(
        chrom = chrom,
        window_width = window_width,
        strand = TRUE
      )
      granges_intersect <- utils.intersect_ranges(granges_dnaseq, granges_binned)
      window_index <- (
        GenomicRanges::findOverlaps( # get bin index for each ribo
          granges_intersect,
          granges_binned,
          select = "first" # only first needed because we intersected
        ) |>
        as.integer()
      )
      # count total coverage per bin
      count <- tapply(
        GenomicRanges::width(granges_intersect) * granges_intersect$count,
        window_index,
        sum
      )
      granges_binned$sample <- S4Vectors::Rle(sample, NROW(granges_binned))
      granges_binned$count <- rep(0, NROW(granges_binned))
      granges_binned$count[as.integer(names(count))] <- count
      
      utils.write_rds(granges_binned, file_out)
    }
  )
}

dnaseq_binned.do_main <- function(
  window_width_list = c(1e3, 1e4, 1e5, 1e6),
  overwrite = FALSE
) {
  purrr::walk(
    window_width_list,
    function(window_width) {
      dnaseq_binned.save(
        window_width = window_width,
        overwrite = overwrite
      )
    }
  )
}
