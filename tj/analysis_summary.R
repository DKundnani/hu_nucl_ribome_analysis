summary.make_sample_table <- function(sample_type, chrom) {
  sample_list <- const.get_sample_list(sample_type)
  chrom_size <- utils.get_chrom_sizes(chrom)

  purrr::map(
    sample_list,
    function(sample) {
      utils.load_ribos_rds(
        sample = sample,
        chrom = chrom
      ) |>
      utils.granges_to_tibble()
    }
  ) |>
  purrr::list_rbind() |>
  dplyr::group_by(sample, chrom, strand) |>
  dplyr::summarize(Ribos = sum(count), .groups = "drop") |>
  dplyr::mutate(
    RPB = sprintf("%.2e", Ribos / chrom_size),
    Ribos = format(Ribos, nsmall = 0, big.mark = ",")
  ) |>
  dplyr::mutate(
    cell = const.get_cell(sample),
    enzyme = const.get_enzyme(sample)
  ) |>
  dplyr::mutate(
    sample = const.get_sample_label_factor(sample, sample_list)
  ) |>
  dplyr::arrange(chrom, strand, sample) |>
  dplyr::select(
    chrom,
    strand,
    Sample = sample,
    Cell = cell,
    Enzyme = enzyme,
    Ribos,
    RPB
  ) |>
  dplyr::nest_by(strand, .keep = TRUE) |>
  purrr::pwalk(
    function(strand, data) {
      data |>
      xtable::xtable(
        caption = paste0(
          "Chrom ", chrom, ".",
          " Chrom size ", as.integer(chrom_size), ".",
          " Strand ", strand, ".",
          " RPB = 'ribos per base'."
        )
      ) |>
      utils.write_xtable(
        fn.summary(
          ext = "tex",
          "sample",
          sample_type,
          chrom,
          strand
        )
      )
    }
  )
}

summary.make_reads_per_window_table <- function(
  sample_type,
  chrom,
  window_width_list
) {
  sample_list <- const.get_sample_list(sample_type)

  purrr::pmap(
    tidyr::expand_grid(
      sample = sample_list,
      window_width = window_width_list
    ),
    function(sample, window_width) {
      utils.load_ribos_binned(
        sample = sample,
        chrom = chrom,
        strand = const.STRAND_BOTH,
        ribo_nuc = const.RIBO_NUC_ALL,
        window_width = window_width,
        normalize_method = "none"
      ) |>
      utils.granges_to_tibble() |>
      dplyr::mutate(
        sample = !!sample,
        window_width = as.integer(!!window_width)
      )
    }
  ) |>
  purrr::list_rbind() |>
  dplyr::nest_by(strand, .keep = TRUE) |>
  purrr::pwalk(
    function(strand, data) {
      data |>
      dplyr::group_by(window_width, sample, chrom, strand) |>
      dplyr::summarize(count_mean = mean(count), .groups = "drop") |>
      dplyr::mutate(window_width = as.integer(window_width)) |>
      tidyr::pivot_wider(
        id_cols = c(sample, chrom, strand),
        names_from = window_width,
        names_sort = TRUE,
        values_from = count_mean
      ) |>
      dplyr::arrange(sample, chrom, strand) |>
      dplyr::rename_with(
        \(x) utils.get_pretty_base_pairs(as.integer(x)),
        .cols = tidyselect::matches("^\\d")
      ) |>
      dplyr::select(!c(chrom, strand)) |>
      dplyr::rename(Sample = sample) |>
      xtable::xtable() |>
      utils.write_xtable(
        fn.summary(
          ext = "tex",
          "reads_per_window",
          sample_type,
          chrom,
          strand
        )
      )
    }
  )
}

summary.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    chrom_list <- "chr1"
    window_width_list <- c(1e5, 1e6)
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ_DRAW
    )
    chrom_list <- const.CHROMS
    window_width_list <- c(1e3, 1e4, 1e5, 1e6)
  }
  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      chrom = chrom_list
    ),
    function(sample_type, chrom) {
      summary.make_sample_table(
        sample_type = sample_type,
        chrom = chrom
      )
      summary.make_reads_per_window_table(
        sample_type = sample_type,
        chrom = chrom,
        window_width_list = window_width_list
      )
    }
  )
}
