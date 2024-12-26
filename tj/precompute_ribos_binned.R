# Get ribo data aggregated into bins across the genome

ribos_binned.save <- function(sample_list, window_width, overwrite = FALSE) {
  purrr::pwalk(
    tidyr::expand_grid(
      sample = sample_list,
      chrom = const.CHROMS
    ),
    function(sample, chrom) {
      file_out <- fn.ribos_binned(
        sample = sample,
        chrom = chrom,
        window_width = window_width
      )
      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }
      ranges_ribo <- utils.load_ribos_rds(sample, chrom)
      if (is.null(ranges_ribo) || (NROW(ranges_ribo) == 0)) {
        utils.log("No ribos for", sample, chrom)
        utils.write_rds(NULL, file_out)
        return(invisible())
      }
      ranges_binned <- utils.get_ranges_binned(
        chrom = chrom,
        window_width = window_width,
        strands = TRUE
      )
      window_index <- (
        GenomicRanges::findOverlaps( # get bin index for each ribo
          ranges_ribo,
          ranges_binned,
          select = "first" # only first needed because ribos are 1-width
        ) |>
        as.integer()
      )

      ranges_binned$sample <- S4Vectors::Rle(sample, NROW(ranges_binned))

      ranges_binned$count <- matrix(
        0,
        nrow = NROW(ranges_binned),
        ncol = length(const.RIBO_NUCS),
        dimnames = list(NULL, const.RIBO_NUCS)
      )
      
      count <- (
        tibble::tibble(
          count = ranges_ribo$count,
          ribo_nuc = ranges_ribo$ribo_nuc,
          window_index = !!window_index
        ) |> # get ribo count and nucleotide for each bin
        dplyr::group_by(ribo_nuc, window_index) |> # group by bin index
        dplyr::summarize(count = sum(count), .groups = "drop") |> # sum ribo counts in each bin
        dplyr::mutate(ribo_nuc = factor(ribo_nuc, levels = const.RIBO_NUCS)) |> # reorder nucleotides
        tidyr::pivot_wider( # pivot to get ribo counts for each nucleotide
          id_cols = window_index,
          names_from = ribo_nuc,
          names_expand = TRUE,
          values_from = count,
          values_fill = 0
        )
      )
      ranges_binned$count[count$window_index, const.RIBO_NUCS] <- (
        as.matrix(count[, const.RIBO_NUCS])
      )
      utils.write_rds(ranges_binned, file_out)
    }
  )
}

ribos_binned.do_main <- function(
  test = FALSE,
  overwrite = FALSE
) {
  if (test) {
    window_width_list <- c(1e5, 1e6)
    sample_list <- c(
      const.SAMPLES_TEST,
      const.SAMPLES_TEST_SHUFFLE,
      const.SAMPLES_TEST_UNIFORM,
      const.SAMPLES_TEST_UNIFORM_CHRALL
    )
  } else {
    window_width_list <- c(1e3, 1e4, 1e5, 1e6)
    sample_list <- c(
      const.SAMPLES,
      const.SAMPLES_SHUFFLE,
      const.SAMPLES_UNIFORM,
      const.SAMPLES_UNIFORM_CHRALL,
      const.SAMPLES_DNASEQ_DRAW
    )
  }
  purrr::walk(
    window_width_list,
    function(window_width) {
      ribos_binned.save(
        sample_list = sample_list,
        window_width = window_width,
        overwrite = overwrite
      )
    }
  )
}
