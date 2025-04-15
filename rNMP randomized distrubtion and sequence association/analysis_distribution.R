# Distribution analysis: Compare the distribution of ribo counts in windows
# between different samples.

distribution.get_tests <- function(x, y, wx, wy, exact, alternative) {
  x <- rep(x, wx)
  y <- rep(y, wy)
  mwu <- wilcox.test(x, y, alternative = alternative, exact = exact)
  ks <- ks.test(x, y, alternative = alternative, exact = exact)
  list(
    mwu_stat = mwu$statistic,
    mwu_p = mwu$p.value,
    ks_stat = ks$statistic,
    ks_p = ks$p.value
  )
}

distribution.get_data <- function(
  sample_type_list, # must have exactly 2 elements: treatment type, control type
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  remove_gaps
) {
  purrr::map(
    sample_type_list,
    function(sample_type) {
      purrr::map(
        const.get_sample_list(sample_type),
        function(sample) {
          utils.load_ribos_binned(
            sample = sample,
            chrom = chrom,
            strand = strand,
            ribo_nuc = ribo_nuc,
            window_width = window_width,
            normalize_method = normalize_method,
            remove_gaps = remove_gaps
          ) |>
          utils.granges_to_tibble() |>
          dplyr::mutate(
            sample = const.get_orig_sample(sample),
            sample_type = sample_type
          ) |>
          dplyr::group_by(
            sample,
            sample_type,
            count
          ) |>
          dplyr::summarize(
            freq = dplyr::n(),
            .groups = "drop"
          )
        }
      ) |>
      purrr::list_rbind()
    }
  ) |>
  purrr::list_rbind()
}

# Do statistical tests comparing 2 types of samples, usually real sample with simulated control.
distribution.make_test_table <- function(
  sample_type_list, # must have exactly 2 elements: treatment type, control type
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  remove_gaps,
  exact, # Use the "exact" setting in the tests.
  alternative # The alternative hypothesis in the tests. Options: "two.sided", "less", "greater"
) {
  sample_list <- (
    const.get_sample_list(sample_type_list[[1]]) |>
    const.get_orig_sample()
  )
  distribution.get_data(
    sample_type_list = sample_type_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width,
    normalize_method = normalize_method,
    remove_gaps = TRUE
  ) |>
  dplyr::nest_by(sample) |>
  purrr::pmap(
    function(data, sample) {
      data_1 <- data |> dplyr::filter(sample_type == sample_type_list[[1]])
      data_2 <- data |> dplyr::filter(sample_type == sample_type_list[[2]])
      test_out <- distribution.get_tests(
        x = data_1[["count"]],
        y = data_2[["count"]],
        wx = data_1[["freq"]],
        wy = data_2[["freq"]],
        exact = exact,
        alternative = alternative
      )
      if (sum(data_1[["freq"]]) != sum(data_2[["freq"]])) {
        stop("The number of windows in the 2 samples are not equal.")
      }
      n_windows <- sum(data_1[["freq"]])
      n_ribos_1 <- sum(data_1[["count"]] * data_1[["freq"]])
      n_ribos_2 <- sum(data_2[["count"]] * data_2[["freq"]])
      mean_1 <- matrixStats::weightedMean(data_1[["count"]], data_1[["freq"]])
      mean_2 <- matrixStats::weightedMean(data_2[["count"]], data_2[["freq"]])
      sd_1 <- matrixStats::weightedSd(data_1[["count"]], data_1[["freq"]])
      sd_2 <- matrixStats::weightedSd(data_2[["count"]], data_2[["freq"]])
      tibble::tibble(
        sample = sample,
        n_windows = n_windows,
        !!paste0("n_ribos_", sample_type_list[[1]]) := n_ribos_1,
        !!paste0("n_ribos_", sample_type_list[[2] ]) := n_ribos_2,
        !!paste0("mean_", sample_type_list[[1]]) := mean_1,
        !!paste0("mean_", sample_type_list[[2]]) := mean_2,
        !!paste0("sd_", sample_type_list[[1]]) := sd_1,
        !!paste0("sd_", sample_type_list[[2]]) := sd_2,
        mwu_stat = test_out[["mwu_stat"]],
        mwu_p = test_out[["mwu_p"]],
        ks_stat = test_out[["ks_stat"]],
        ks_p = test_out[["ks_p"]]
      )
    }
  ) |>
  purrr::list_rbind() |>
  utils.write_csv(
    fn.distribution_table(
      prefix = "test",
      sample_type_list = sample_type_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_width = window_width,
      normalize_method = normalize_method,
      remove_gaps = remove_gaps,
      if (exact) "exact" else "notExact",
      stringr::str_replace(alternative, "\\.", ""),
      ext = "csv"
    )
  )
}

# Density plot of ribo counts in windows comparing 2 types of samples.
distribution.plot_density <- function(
  sample_type_list, # must have exactly 2 elements: treatment type, control type
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  remove_gaps,
  # Percentile cutoff for the x-axis.
  # The x-axis will be truncated at this percentile of density for the first sample type.
  p_cutoff,
  base_font_size = 11,
  sample_type_name_list = NULL, # The labels to use for the sample types.
  remove_legend = FALSE, # Whether to remove the legend.
  remove_axis_labels = FALSE, # Whether to remove the axis labels.
  n_ticks_y = NULL # Number of ticks for the y-axis
) {
  if (is.null(sample_type_name_list)) {
    sample_type_name_list <- sample_type_list
  }

  distribution.get_data(
    sample_type_list = sample_type_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width,
    normalize_method = normalize_method,
    remove_gaps = remove_gaps
  ) |>
  dplyr::mutate(
    sample_name = factor(
      sample_type,
      levels = sample_type_list,
      labels = sample_type_name_list
    )
  ) |>
  dplyr::nest_by(sample) |>
  purrr::pwalk(
    function(data, sample) {
      cutoff <- ( # x-axis cutoff
        data |>
        dplyr::filter(sample_type == sample_type_list[[1]]) |>
        x => utils.get_quantiles(x[["count"]], x[["freq"]], p_cutoff)
      )

      # Maks sure both samples have all x-axis values
      max_val <- data |> dplyr::pull(count) |> max()
      data <- (
        data |>
        dplyr::right_join(
          tidyr::expand_grid(
            count = seq(0, max_val),
            sample_type = sample_type_list
          ),
          by = c("count", "sample_type")
        ) |>
        tidyr::replace_na(list(freq = 0))
      )

      (
        data |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = count,
            y = freq,
            fill = sample_name
          )
        ) +
        ggplot2::geom_col(width = 1, position = "dodge") +
        ggplot2::coord_cartesian(xlim = c(0, cutoff)) +
        ggplot2::scale_y_continuous(n.breaks = n_ticks_y) +
        ggplot2::theme_classic(base_size = base_font_size) +
        ggplot2::xlab("Ribo count in window") +
        ggplot2::ylab("Number of windows") +
        ggplot2::labs(
          title = stringr::str_c(
            "Sample: ", sample, "\n",
            sample_type_list[[1]], " vs. ", sample_type_list[[2]], " distributions\n",
            "Chromsome: ", chrom, "\n",
            "Strand: ", strand, "\n",
            "Ribo nucs: ", ribo_nuc, "\n",
            "Window width: ", utils.get_pretty_base_pairs(window_width), "\n",
            "Normalize method: ", normalize_method, "\n",
            "Percentile cutoff: ", (p_cutoff * 100), "\n",
            "Y-axis shows the number of windows\n",
            "X-axis shows the ribo count in the window"
          )
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          legend.position = if (remove_legend) "none" else "top",
          axis.title = if (remove_axis_labels) {
            ggplot2::element_blank()
          } else {
            ggplot2::element_text(size = base_font_size)
          },
          legend.title = ggplot2::element_blank(),
          # Made right margins 2 cm to make space for the tick labels
          plot.margin = ggplot2::margin(t = 1, r = 2, b = 1, l = 1, "cm")
        )
      ) |>
      utils.write_ggplot(
        fn.distribution_plot(
          prefix = "density",
          sample_type_list = sample_type_list,
          chrom = chrom,
          strand = strand,
          ribo_nuc = ribo_nuc,
          window_width = window_width,
          normalize_method = normalize_method,
          remove_gaps = remove_gaps,
          p_cutoff = round(p_cutoff * 100),
          sample = sample
        ),
        width = 16,
        height = 4,
        limitsize = FALSE
      )
    }
  )
}

distribution.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list_list <- list(
      c(const.SAMPLE_TEST, const.SAMPLE_TEST_SHUFFLE)
    )
    chrom_list <- "chr1"
    window_width_list <- 1e5
    strand <- const.STRAND_BOTH
    ribo_nuc <- const.RIBO_NUC_ALL
  } else {
    sample_type_list_list <- list(
      c(const.SAMPLE_NORMAL, const.SAMPLE_SHUFFLE),
      c(const.SAMPLE_NORMAL, const.SAMPLE_UNIFORM),
      c(const.SAMPLE_NORMAL, const.SAMPLE_UNIFORM_CHRALL),
      c(const.SAMPLE_NORMAL, const.SAMPLE_DNASEQ_DRAW)
    )
    chrom_list <- const.CHROMS
    window_width_list <- c(1e3, 1e4, 1e5, 1e6)
    strand <- const.STRAND_BOTH
    ribo_nuc <- const.RIBO_NUC_ALL
    p_cutoff <- 0.99
  }

  for (chrom in chrom_list) {
    for (sample_type_list in sample_type_list_list) {
      for (window_width in window_width_list) {
        distribution.plot_density(
          sample_type_list = sample_type_list,
          chrom = chrom,
          strand = strand,
          ribo_nuc = ribo_nuc,
          window_width = window_width,
          normalize_method = "none",
          remove_gaps = TRUE,
          p_cutoff = p_cutoff,
          base_font_size = settings.BASE_FONT_SIZE
        )
      }
    }
  }
}

