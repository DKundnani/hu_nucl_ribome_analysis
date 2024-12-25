overlap_clusters.get_data <- function(
  sample,
  chrom,
  window_width_list
) {
  ribos <- (
    utils.load_ribos_rds(sample = sample, chrom = chrom) |>
    utils.granges_to_tibble()
  )

  tidyr::expand_grid(
    strand = const.STRANDS,
    window_width = window_width_list
  ) |>
  purrr::pmap(
    function(strand, window_width) {
      ribos |> dplyr::filter(strand == !!strand) |>
      x => utils.reduce_intervals(
        x[["start"]] - window_width,
        x[["start"]] + window_width
      ) |>
      dplyr::mutate(
        start = start + window_width,
        end = end - window_width
      ) |>
      dplyr::mutate(width = end - start + 1) |>
      dplyr::summarize(
        sample = !!sample,
        chrom = !!chrom,
        strand = !!strand,
        window_width = !!window_width,
        cluster_n = dplyr::n(),
        cluster_n_2 = sum(width >= 2),
        cluster_width_mean_2 = mean(width[width >= 2]),
        cluster_width_median_2 = median(width[width >= 2]),
        cluster_width_max_2 = max(width[width >= 2])
      )
    }
  ) |>
  purrr::list_rbind()
}

overlap_clusters.plot <- function(
  sample_type,
  chrom,
  window_width_list,
  base_font_size = 11
) {
  sample_list = const.get_sample_list(sample_type)

  data <- (
    purrr::map(
      sample_list,
      function(sample) {
        overlap_clusters.get_data(
          sample = sample,
          chrom = chrom,
          window_width_list = window_width_list
        )
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::mutate(window_width = factor(window_width, window_width_list))
  )

  data |>
  dplyr::nest_by(strand, .key = "data_strand") |>
  tidyr::expand_grid(
    col_name = stringr::str_subset(colnames(data), "cluster_")
  ) |>
  purrr::pwalk(
    function(strand, data_strand, col_name) {
      (
        data_strand |>
        dplyr::mutate(
          sample = const.get_sample_label_factor(sample, sample_list)
        ) |>
        ggplot2::ggplot(
          mapping = ggplot2::aes(
            x = window_width,
            y = get(col_name)
          )
        ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = stringr::str_c(
            "Clusters summary expanding by window_width and merging\n",
            "Summary statistic = ", col_name, "\n",
            "Sample type = ", sample_type, "\n",
            "Chromosome = ", chrom, "\n",
            "Strand = ", strand, "\n",
            "X-axis show window width for merging clusters\n",
            "Y-axis shows summary statistic ", col_name, " for merged clusters\n",
            "Each panel shows a different sample"
          )
        ) +
        ggplot2::xlab("Window size (bp)") +
        ggplot2::ylab(col_name) +
        ggplot2::theme_classic(base_size = base_font_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          strip.text = ggplot2::element_text(size = 7, angle = 0),
          strip.background = ggplot2::element_rect(color = "black"),
          legend.position = "top",
          legend.text = ggplot2::element_text(),
          legend.title = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(),
          axis.text.x = ggplot2::element_text(angle = 90)
        ) +
        ggplot2::facet_wrap(
          dplyr::vars(sample),
          scales = "free_y"
        )
      ) |>
      utils.write_ggplot(
        fn.overlap_clusters_plot(
          sample_type = sample_type,
          chrom = chrom,
          strand = strand,
          col_name = col_name
        ),
        width = (
          3 *
          ceiling(sqrt(length(sample_list))) *
          max(length(window_width_list) / 5, 1)
        ),
        height = 2 * ceiling(sqrt(length(sample_list))),
        limitsize = FALSE
      )
    }
  )
}

overlap_clusters.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    chrom_list <- "chr1"
    window_width_list <- c(0, 1, 10, 100)
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ_DRAW
    )
    chrom_list <- const.CHROMS
    window_width_list <- c(0, 1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 75, 100)
  }

  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      chrom = chrom_list
    ),
    function(sample_type, chrom) {
      overlap_clusters.plot(
        sample_type = sample_type,
        chrom = chrom,
        window_width_list = window_width_list,
        base_font_size = settings.BASE_FONT_SIZE
      )
    }
  )
}
