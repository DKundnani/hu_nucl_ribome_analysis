between_ribo_dist.get_data <- function(
  sample_list,
  chrom,
  diff_breaks,
  prob_type # choices: "pdf", "cdf"
) {
  diff_labels <- purrr::map2_chr(
    diff_breaks[-length(diff_breaks)],
    diff_breaks[-1],
    function(start, end) {
      if (is.infinite(end)) {
        paste0(start, "+")
      } else if (start == end) {
        as.character(start)
      } else {
        paste0(start, " to ", end - 1)
      }
    }
  )

  bin_diffs <- function(pos) {
    pos |>
    sort() |>
    diff() |>
    cut(breaks = diff_breaks, labels = diff_labels, right = FALSE) |>
    as.character() |>
    x => tibble::tibble(diff = x) |>
    dplyr::group_by(diff) |>
    dplyr::summarize(freq = dplyr::n()) |>
    dplyr::right_join(tibble::tibble(diff = diff_labels), by = "diff") |>
    dplyr::mutate(freq = tidyr::replace_na(freq, 0)) |>
    dplyr::mutate(diff = factor(diff, diff_labels))
  }

  purrr::map(
    sample_list,
    function(sample) {
      utils.load_ribos_rds(sample = sample, chrom = chrom) |>
      utils.granges_to_tibble() |>
      dplyr::nest_by(strand) |>
      purrr::pmap(
        function(strand, data) {
          data <- (
            data |>
            dplyr::pull(start) |>
            bin_diffs()
          )
          if (prob_type == "cdf") {
            data <- (
              data |>
              dplyr::arrange(diff) |>
              dplyr::mutate(freq = cumsum(freq)) |>
              dplyr::mutate(
                diff = forcats::fct_relabel(
                  diff,
                  function(x) stringr::str_c("<", diff_breaks[-1])
                )
              )
            )
          }

          data |>
          dplyr::mutate(
            sample = !!sample,
            chrom = !!chrom,
            strand = !!strand,
            .before = 1
          )
        }
      ) |> purrr::list_rbind()
    }
  ) |>
  purrr::list_rbind()
}

between_ribo_dist.plot <- function(
  sample_type_list,
  chrom,
  diff_breaks,
  prob_type, # choices: "pdf", "cdf"
  base_font_size = 11
) {
  num_samples <- (
    purrr::map_int(
      sample_type_list,
      function(sample_type) {
        length(const.get_sample_list(sample_type))
      }
    ) |>
    max()
  )
  data <- (
    purrr::map(
      sample_type_list,
      function(sample_type) {
        between_ribo_dist.get_data(
          sample_list = const.get_sample_list(sample_type),
          chrom = chrom,
          diff_breaks = diff_breaks,
          prob_type = prob_type
        ) |>
        dplyr::mutate(
          sample_type = !!sample_type,
          sample = const.get_orig_sample(sample),
          .before = 1
        )
      }
    ) |>
    purrr::list_rbind()
  )

  data |>
  dplyr::nest_by(strand, .key = "data_strand") |>
  purrr::pwalk(
    function(strand, data_strand) {
      (
        data_strand |>
        dplyr::mutate(
          sample = const.get_sample_label_factor(
            sample,
            const.get_sample_list(const.SAMPLE_NORMAL)
          )
        ) |>
        ggplot2::ggplot(
          mapping = ggplot2::aes(x = diff, y = freq, fill = sample_type)
        ) +
        ggplot2::geom_col(position = ggplot2::position_dodge()) +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
          title = paste0(
            "Between ribo distance distribution\n",
            "Sample types = ", paste0(sample_type_list, collapse = ", "), "\n",
            "Chromosome = ", chrom, "\n",
            "Strand = ", strand, "\n",
            "Probability type = ", prob_type, "\n",
            "Columns show means and error bars show min/max\n",
            "Each panel shows a different sample"
          )
        ) +
        ggplot2::xlab("Difference (bp)") +
        ggplot2::ylab(
          if (prob_type == "pdf") {
            "Frequency (count)"
          } else {
            "Frequency (cumulative count)"
          }
        ) +
        ggplot2::theme_classic(base_size = base_font_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          strip.text = ggplot2::element_text(size = 7, angle = 0),
          strip.background = ggplot2::element_rect(color = "black"),
          legend.position = "top",
          legend.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 90)
        ) +
        ggplot2::facet_wrap(dplyr::vars(sample), scales = "free_y")
      ) |>
      utils.write_ggplot(
        fn.between_ribo_dist_plot(
          sample_type = paste0(sample_type_list, collapse = "_"),
          chrom = chrom,
          strand = strand,
          prob_type = prob_type
        ),
        width = max(ceiling(2 * sqrt(num_samples) * length(diff_breaks) / 5), 4),
        height = max(ceiling(2 * sqrt(num_samples)), 4),
        limitsize = FALSE
      )
    }
  )
}

between_ribo_dist.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list_list <- list(
      c(const.SAMPLE_TEST, const.SAMPLE_TEST_SHUFFLE)
    )
    chrom_list <- "chr1"
  } else {
    sample_type_list_list <- list(
      c(const.SAMPLE_NORMAL, const.SAMPLE_SHUFFLE),
      c(const.SAMPLE_NORMAL, const.SAMPLE_UNIFORM),
      c(const.SAMPLE_NORMAL, const.SAMPLE_UNIFORM_CHRALL),
      c(const.SAMPLE_NORMAL, const.SAMPLE_DNASEQ_DRAW)
    )
    chrom_list <- const.CHROMS
  }
  prob_type_list <- c("pdf", "cdf")
  diff_breaks <- c(1, 10, 25, 50, 100, 500, 1000, Inf)

  purrr::pwalk(
    tidyr::expand_grid(
      sample_type_list = sample_type_list_list,
      chrom = chrom_list,
      prob_type = prob_type_list
    ),
    function(sample_type_list, chrom, prob_type) {
      between_ribo_dist.plot(
        sample_type_list = sample_type_list,
        chrom = chrom,
        diff_breaks = diff_breaks,
        prob_type = prob_type,
        base_font_size = settings.BASE_FONT_SIZE
      )
    }
  )
}
