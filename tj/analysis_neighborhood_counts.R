neighborhood_counts.get_data <- function(
  sample_list,
  chrom,
  neighborhood_size,
  count_type
) {
  purrr::map(
    sample_list,
    function(sample) {
      utils.get_neighborhood_counts(
        sample = sample,
        chrom = chrom,
        neighborhood_size = neighborhood_size,
        count_type = count_type,
        only_orig = TRUE
      ) |>
      dplyr::mutate(sample = !!sample)
    }
  ) |>
  purrr::list_rbind() |>
  dplyr::group_by(sample, chrom, strand, neighborhood_size, count_nb) |> # Get the empirical distribution of counts
  dplyr::summarize(freq = dplyr::n(), .groups = "drop") |>
  dplyr::rename(count = count_nb)
}

neighborhood_counts.plot <- function(
  sample_type,
  chrom,
  neighborhood_size,
  count_type,
  max_count = NULL,
  base_font_size = 11
) {
  sample_list <- const.get_sample_list(sample_type)

  neighborhood_counts.get_data(
    sample_list = sample_list,
    chrom = chrom,
    neighborhood_size = neighborhood_size,
    count_type = count_type
  ) |>
  dplyr::nest_by(strand) |>
  purrr::pwalk(
    function(strand, data) {
      (
        data |>
        dplyr::mutate(
          sample = const.get_sample_label_factor(sample, sample_list)
        ) |>
        ggplot2::ggplot() +
        ggplot2::geom_col(
          mapping = ggplot2::aes(x = count, y = freq)
        ) +
        ggplot2::scale_y_continuous(trans = "log10") +
        (
          if (!is.null(max_count)) {
            ggplot2::expand_limits(x = max_count)
          } else {
            ggplot2::geom_blank()
          }
        ) +
        (
          if (!is.null(max_count)) {
            ggplot2::xlim(c(0, max_count))
          } else {
            ggplot2::geom_blank()
          }
        ) +
        ggplot2::labs(
          title = paste0(
            "Neighborhood ribo count histogram\n",
            "Chromosome = ", chrom, "\n",
            "Strand = ", strand, "\n",
            "Neighborhood size = ", utils.get_pretty_base_pairs(as.integer(neighborhood_size)), "\n",
            "Count type = ", count_type, "\n",
            "Sample type = ", sample_type, "\n",
            "X-axis shows the number of ribos in a neighborhood\n",
            "Y-axis shows the number of neighborhoods with that many ribos"
          )
        ) +
        ggplot2::xlab(paste0("Neighborhood ribo count (", count_type, ")")) +
        ggplot2::ylab("Number of distinct genome loci (log scale)") +
        ggplot2::facet_wrap(
          dplyr::vars(sample),
          scales = "free"
        ) +
        ggplot2::theme_classic(base_size = base_font_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          strip.text = ggplot2::element_text(size = 7, angle = 0),
          strip.background = element_rect(color = "black"),
          legend.position = "none",
          legend.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 90)
        )
      ) |>
      utils.write_ggplot(
        fn.neighborhood_counts_plot(
          sample_type = sample_type,
          chrom = chrom,
          strand = strand,
          neighborhood_size = neighborhood_size,
          count_type = count_type
        ),
        width = ceiling(sqrt(length(sample_list))) * 3,
        height = ceiling(sqrt(length(sample_list))) * 3
      )
    }
  )
}

neighborhood_counts.table <- function(
  sample_type,
  chrom,
  neighborhood_size,
  count_type,
  count_breaks
) {
  sample_list <- const.get_sample_list(sample_type)
  neighborhood_counts.get_data(
    sample_list = sample_list,
    chrom = chrom,
    neighborhood_size = neighborhood_size,
    count_type = count_type
  ) |>
  dplyr::mutate(
    sample = purrr::map_chr(as.character(sample), const.get_label),
    count = cut(
      count,
      breaks = count_breaks,
      labels = purrr::map_chr(
        seq_len(length(count_breaks) - 1),
        function(i) {
          if (count_breaks[[i]] == (count_breaks[[i + 1]] - 1)) {
            as.character(count_breaks[[i]])
          } else {
            paste0(count_breaks[[i]], " to ", count_breaks[[i + 1]] - 1)
          }
        }
      ),
      right = FALSE
    )
  ) |>
  dplyr::group_by(sample, chrom, strand, count) |>
  dplyr::summarize(freq = as.integer(sum(freq)), .groups = "drop") |>
  dplyr::nest_by(strand, .keep = TRUE) |>
  purrr::pwalk(
    function(strand, data) {
      data |>
      tidyr::pivot_wider(
        id_cols = sample,
        names_from = count,
        values_from = freq,
        values_fill = 0
      ) |>
      xtable::xtable() |>
      utils.write_xtable(
        fn.neighborhood_counts_table(
          sample_type = sample_type,
          chrom = chrom,
          strand = strand,
          neighborhood_size = neighborhood_size,
          count_type = count_type
        )
      )
    }
  )
}

neighborhood_counts.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    chrom_list <- "chr1"
    neighborhood_size_list <- 500L
    count_type_list <- c("multiple", "single")
    output_type_list <- c("plot", "table")
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ
    )
    chrom_list <- const.CHROMS
    neighborhood_size_list <- c(0L, 10L, 50L, 100L, 500L, 1000L)
    count_type_list <- c("multiple", "single")
    output_type_list <- c("plot", "table")
  }

  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      chrom = chrom_list,
      neighborhood_size = neighborhood_size_list,
      count_type = count_type_list,
      output_type = output_type_list
    ),
    function(
      sample_type,
      chrom,
      neighborhood_size,
      count_type,
      output_type
    ) {
      if (output_type == "plot") {
        neighborhood_counts.plot(
          sample_type = sample_type,
          chrom = chrom,
          neighborhood_size = neighborhood_size,
          count_type = count_type,
          max_count = NULL,
          base_font_size = settings.BASE_FONT_SIZE
        )
      } else if (output_type == "table") {
        neighborhood_counts.table(
          sample_type = sample_type,
          chrom = chrom,
          neighborhood_size = neighborhood_size,
          count_type = count_type,
          count_breaks = c(0L, 1L, 2L, 3L, 4L, 5L, 10L, 20L, 40L, 100L, Inf)
        )
      } else {
        stop(paste0("Unknown output type: ", output))
      }
    }
  )
}
