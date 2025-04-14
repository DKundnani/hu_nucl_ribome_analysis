library(transport) # for wasserstein1d
library(dtw) # for dtw

line_graph.OUTLIER_INDICATOR_FUNCS <- list(
  "p_value" = function(x, p_value) {
    ifelse(x >= quantile(x, 1 - p_value), 1, 0)
  },
  "sd_thresh" = function(x, sd_thresh) {
    ifelse((x - mean(x)) >= sd_thresh * sd(x), 1, 0)
  }
)

line_graph.format_args <- function(args) {
  paste(names(args), "=", args, collapse = ", ")
}

line_graph.get_outlier_indicator <- function(x, outlier_type, outlier_args) {
  do.call(
    line_graph.OUTLIER_INDICATOR_FUNCS[[outlier_type]],
    args = c(list(x = x), outlier_args)
  )
}

line_graph.get_outlier_indicator_data <- function(
  data,
  outlier_type,
  outlier_args
) {
  data |>
  dplyr::group_by(sample, strand) |>
  dplyr::mutate(
    outlier = line_graph.get_outlier_indicator(
      x = count,
      outlier_type = outlier_type,
      outlier_args = outlier_args
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(strand, start) |>
  dplyr::summarize(n = sum(outlier), .groups = "drop")
}

line_graph.SCORE_TYPES <- list(
  sim_jacc_01 = list(
    func = function(x, y) {
      x <- ifelse(line_graph.get_outlier_indicator(x, 0.01), 1, 0)
      y <- ifelse(line_graph.get_outlier_indicator(y, 0.01), 1, 0)
      sum(x & y) / sum(x | y)
    },
    label = "Jaccard similarity of 1% outlier windows",
    limits = c(0, 1),
    breaks = c(0, 1),
    colors = c("white", "black")
  ),
  sim_jacc_05 = list(
    func = function(x, y) {
      x <- ifelse(line_graph.get_outlier_indicator(x, 0.05), 1, 0)
      y <- ifelse(line_graph.get_outlier_indicator(y, 0.05), 1, 0)
      sum(x & y) / sum(x | y)
    },
    label = "Jaccard similarity of 5% outlier windows",
    limits = c(0, 1),
    breaks = c(0, 1),
    colors = c("white", "black")
  ),
  sim_jacc_10 = list(
    func = function(x, y) {
      x <- ifelse(line_graph.get_outlier_indicator(x, 0.1), 1, 0)
      y <- ifelse(line_graph.get_outlier_indicator(y, 0.1), 1, 0)
      sum(x & y) / sum(x | y)
    },
    label = "Jaccard similarity of 10% outlier windows",
    limits = c(0, 1),
    breaks = c(0, 1),
    colors = c("white", "black")
  ),
  corr_pearson = list(
    func = function(x, y) cor(x, y, method = "pearson"),
    label = "Pearson correlation",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    colors = c("red", "white", "blue")
  ),
  corr_spearman = list(
    func = function(x, y) cor(x, y, method = "spearman"),
    label = "Spearman correlation",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    colors = c("red", "white", "blue")
  ),
  corr_kendall = list(
    func = function(x, y) cor(x, y, method = "kendall"),
    label = "Kendall correlation",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    colors = c("red", "white", "blue")
  ),
  dist_l1 = list(
    func = function(x, y) {
      d <- function(u) sum(abs(u))
      x <- x / d(x)
      y <- y / d(y)
      d(x - y)
    },
    label = "L^1 distance",
    limits = c(0, 2),
    breaks = c(0, 2),
    colors = c("white", "black")
  ),
  dist_l2 = list(
    func = function(x, y) {
      d <- function(u) sqrt(sum(abs(u)^2))
      x <- x / d(x)
      y <- y / d(y)
      d(x - y)
    },
    label = "L^2 distance",
    limits = c(0, 2),
    breaks = c(0, 2),
    colors = c("white", "black")
  ),
  dist_l3 = list(
    func = function(x, y) {
      d <- function(u) sum(abs(u)^3)^(1/3)
      x <- x / d(x)
      y <- y / d(y)
      d(x - y)
    },
    label = "L^3 distance",
    limits = c(0, 2),
    breaks = c(0, 2),
    colors = c("white", "black")
  ),
  dist_max = list(
    func = function(x, y) {
      d <- function(u) max(u)
      x <- x / d(x)
      y <- y / d(y)
      d(x - y)
    },
    label = "Max distance",
    limits = c(0, 2),
    breaks = c(0, 2),
    colors = c("white", "black")
  ),
  dist_wasser = list(
    func = function(x, y) {
      transport::wasserstein1d(
        a = seq_along(x),
        b = seq_along(x),
        p = 1,
        wa = x,
        wb = y
      )
    },
    label = "Wasserstein^1 distance",
    limits = NULL,
    breaks = NULL,
    colors = NULL
  ),
  dist_time_warp = list(
    func = function(x, y) {
      x <- x / sqrt(sum(x^2))
      y <- y / sqrt(sum(y^2))
      tryCatch(
        dtw::dtw(x, y)[["distance"]],
        error = function(e) NA
      )
    },
    label = "Dynamic-time-warp distance",
    limits = NULL,
    breaks = NULL,
    colors = NULL
  )
)

line_graph.get_data <- function(
  sample_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method
) {
  purrr::map(
    sample_list,
    function(sample) {
      utils.load_ribos_binned(
        sample = sample,
        chrom = chrom,
        strand = strand,
        ribo_nuc = ribo_nuc,
        window_width = window_width,
        normalize_method = normalize_method
      ) |>
      utils.granges_to_tibble() |>
      dplyr::mutate(sample = sample)
    }
  ) |>
  purrr::list_rbind()
}

line_graph.plot_line <- function(
  sample_type,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  base_font_size = 11
) {
  if (chrom == const.CHROM_ALL) {
    stop("Cannot plot all chromosomes on a line graph.")
  }
  if (strand == const.STRAND_BOTH) {
    stop("Cannot plot both strands on a line graph.")
  }

  sample_list <- const.get_sample_list(sample_type)
  tick_info <- utils.get_tick_info(utils.get_chrom_sizes(chrom))

  (
    line_graph.get_data(
      sample_list = sample_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_width = window_width,
      normalize_method = normalize_method
    ) |>
    dplyr::mutate(
      sample = const.get_sample_label_factor(sample, sample_list)
    ) |>
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = start,
        y = count,
        group = sample
      )
    ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_x_continuous(
      breaks = tick_info[["breaks"]],
      labels = tick_info[["labels"]]
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::ylab("Number of reads") +
    ggplot2::xlab("Chromosome position") +
    ggplot2::labs(
      title = paste0(
        "Binned count line graph\n",
        "Sample type = ", sample_type, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nucs = ", ribo_nuc, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "Normalize method = ", normalize_method, "\n",
        "X-axis shows the position on the chromosome\n",
        "Y-axis shows the read count in the window\n",
        "Each panel is a different sample"
      )
    ) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      strip.placement = "top",
      strip.text.y = ggplot2::element_text(size = 7, angle = 0),
      strip.background = ggplot2::element_rect(color = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90)
    ) +
    ggplot2::facet_grid(
      rows = dplyr::vars(sample),
      scales = "free_y"
    )
  ) |>
  utils.write_ggplot(
    fn.line_graph_plot(
      sample_type = sample_type,
      window_width = window_width,
      dir = file.path("line_graph", sample_type),
      prefix = "line_graph",
      chrom,
      strand,
      ribo_nuc,
      normalize_method
    ),
    width = 12 + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = max(2 * length(sample_list), 4),
    limitsize = FALSE
  )
}

line_graph.plot_outlier <- function(
  sample_type,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  outlier_type,
  outlier_args,
  base_font_size = 11
) {
  if (chrom == const.CHROM_ALL) {
    stop("Cannot plot outliers on all chromosomes.")
  }
  if (strand == const.STRAND_BOTH) {
    stop("Cannot plot outliers on both strands.")
  }

  sample_list <- const.get_sample_list(sample_type)
  chrom_size <- utils.get_chrom_sizes(chrom)
  tick_info <- utils.get_tick_info(chrom_size)

  (
    line_graph.get_data(
      sample_list = sample_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_width = window_width,
      normalize_method = normalize_method
    ) |>
    tidyr::drop_na(count) |>
    dplyr::mutate(
      sample = const.get_sample_label_factor(sample, sample_list)
    ) |>
    dplyr::group_by(sample) |>
    dplyr::mutate(
      outlier = line_graph.get_outlier_indicator(
        x = count,
        outlier_type = outlier_type,
        outlier_args = outlier_args
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(outlier != 0) |>
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = start,
        xend = start + window_width,
        y = sample,
        yend = sample,
        group = sample
      )
    ) +
    ggplot2::geom_segment(linewidth = 10) +
    ggplot2::scale_x_continuous(
      breaks = tick_info[["breaks"]],
      labels = tick_info[["labels"]]
    ) +
    ggplot2::ylab("Sample") +
    ggplot2::xlab("Chromosome position") +
    ggplot2::labs(
      title = paste0(
        "Outlier window segments\n",
        "Sample type = ", sample_type, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nucs = ", ribo_nuc, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "Normalize method = ", normalize_method, "\n",
        "Outlier type = ", outlier_type, "\n",
        "Outlier args = ", line_graph.format_args(outlier_args), "\n",
        "X-axis shows the position on the chromosome\n",
        "Y-axis shows the sample"
      )
    ) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      strip.placement = "top",
      strip.text.x = ggplot2::element_text(angle = 0),
      strip.text.y = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_rect(color = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90)
    )
  ) |>
  utils.write_ggplot(
    fn.line_graph_plot(
      sample_type = sample_type,
      window_width = window_width,
      dir = file.path("outlier_windows", sample_type),
      prefix = "outlier_windows",
      as.character(chrom),
      as.character(strand),
      as.character(ribo_nuc),
      normalize_method,
      outlier_type,
      paste(outlier_args, collapse = "_")
    ),
    width = 12 + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = length(sample_list),
    limitsize = FALSE
  )
}

line_graph.plot_outlier_indicators <- function(
  sample_type,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  outlier_type,
  outlier_args,
  base_font_size = 11
) {
  if (chrom == const.CHROM_ALL) {
    stop("Cannot plot outlier indicators for all chromosomes.")
  }
  if (strand == const.STRAND_BOTH) {
    stop("Cannot plot outlier indicators on both strands.")
  }

  sample_list <- const.get_sample_list(sample_type)
  chrom_size <- utils.get_chrom_sizes(chrom)
  tick_info <- utils.get_tick_info(chrom_size)
  
  (
    line_graph.get_data(
      sample_list = sample_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_width = window_width,
      normalize_method = normalize_method
    ) |>
    tidyr::drop_na(count) |>
    line_graph.get_outlier_indicator_data(
      outlier_type = outlier_type,
      outlier_args = outlier_args
    ) |>
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = start,
        y = n
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(
      breaks = tick_info[["breaks"]],
      labels = tick_info[["labels"]]
    ) +
    ggplot2::ylab("Number of samples") +
    ggplot2::xlab("Chromosome position") +
    ggplot2::labs(
      title = paste0(
        "Outlier indicator sum\n",
        "Sample type = ", sample_type, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nuc = ", ribo_nuc, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "Normalize method = ", normalize_method, "\n",
        "Outlier type = ", outlier_type, "\n",
        "Outlier args = ", line_graph.format_args(outlier_args), "\n",
        "X-axis shows the position on the chromosome\n",
        "Y-axis shows the number of samples with an outlier window"
      )
    ) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      strip.text.x = ggplot2::element_text(angle = 0),
      strip.text.y = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_rect(color = "black"),
      axis.text.x = ggplot2::element_text(angle = 90)
    )
  ) |>
  utils.write_ggplot(
    fn.line_graph_plot(
      sample_type = sample_type,
      window_width = window_width,
      dir = file.path("outlier_indicators", sample_type),
      prefix = "outlier_indicators",
      chrom,
      strand,
      ribo_nuc,
      normalize_method,
      outlier_type,
      paste(outlier_args, collapse = "_")
    ),
    width = 12 + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = 3,
    limitsize = FALSE
  )
}

line_graph.plot_matrix <- function(
  sample_type,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  score_type,
  base_font_size
) {
  if (!(score_type %in% names(line_graph.SCORE_TYPES))) {
    stop(paste0("Unknown score type: ", score_type))
  }
  sample_list <- const.get_sample_list(sample_type)
  count_data <- (
    line_graph.get_data(
      sample_list = sample_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_width = window_width,
      normalize_method = normalize_method
    ) |>
    tidyr::drop_na(count) |>
    dplyr::arrange(sample, chrom, strand, start) |>
    dplyr::nest_by(sample)
  )

  matrix_data <- (
    purrr::map(
      combn(nrow(count_data), 2, simplify = FALSE),
      function(index) {
        sample_1 <- count_data[[index[[1]], "sample"]]
        sample_2 <- count_data[[index[[2]], "sample"]]
        count_1 <- count_data[["data"]][[index[[1]]]][["count"]]
        count_2 <- count_data[["data"]][[index[[2]]]][["count"]]
        utils.log(strand, " ", as.character(sample_1), " ", as.character(sample_2))
        score <- line_graph.SCORE_TYPES[[score_type]][["func"]](count_1, count_2)
        tibble::tibble(
          strand = !!strand,
          sample_1 = c(!!sample_1, !!sample_2),
          sample_2 = c(!!sample_2, !!sample_1),
          score = c(!!score, !!score)
        )
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::bind_rows(
      purrr::map(
        seq_len(nrow(count_data)),
        function(index) {
          sample <- count_data[["sample"]][[index]]
          count <- count_data[["data"]][[index]][["count"]]
          score <- line_graph.SCORE_TYPES[[score_type]][["func"]](count, count)
          tibble::tibble(
            strand = !!strand,
            sample_1 = sample,
            sample_2 = sample,
            score = !!score
          )
        }
      ) |>
      purrr::list_rbind()
    )
  )

  (
    matrix_data |>
    dplyr::mutate(
      score_label = format(round(score, 2), digits = 2),
      sample_1 = const.get_sample_label_factor(sample_1, sample_list),
      sample_2 = const.get_sample_label_factor(sample_2, sample_list)
    ) |>
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = sample_1,
        y = sample_2,
        fill = score,
        label = score_label
      )
    ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text() +
    ggplot2::labs(
      title = paste0(
        "Line graph pairwise score matrix\n",
        "Sample type = ", sample_type, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nucs = ", ribo_nuc, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "Normalize method = ", normalize_method, "\n",
        "X-axis/y-axis show the samples being compared\n",
        "Color shows the ",
        line_graph.SCORE_TYPES[[score_type]][["label"]], " score\n",
        "Light pink is N/A\n"
      )
    ) +
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      strip.text = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_rect(color = "black"),
      legend.position = "right",
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      axis.text.y = ggplot2::element_text(hjust = 1),
      aspect.ratio = 1
    ) +
    ggplot2::scale_fill_gradientn(
      name = line_graph.SCORE_TYPES[[score_type]][["label"]],
      limits = line_graph.SCORE_TYPES[[score_type]][["limits"]],
      breaks = line_graph.SCORE_TYPES[[score_type]][["breaks"]],
      colors = line_graph.SCORE_TYPES[[score_type]][["colors"]],
      na.value = "#fff0f0",
      oob = scales::squish
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barheight = ggplot2::unit(3, "in"),
        barwidth = ggplot2::unit(0.5, "in"),
        ticks.colour = "black",
        ticks.linewidth = 1,
        frame.colour = "black",
        frame.linewidth = 1
      )
    )
  ) |>
  utils.write_ggplot(
    fn.line_graph_plot(
      sample_type = sample_type,
      window_width = window_width,
      dir = file.path("matrix", score_type, sample_type),
      prefix = "matrix",
      score_type,
      chrom,
      strand,
      ribo_nuc,
      normalize_method
    ),
    width = length(sample_list) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = length(sample_list) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    limitsize = FALSE
  )
}

line_graph.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    chrom_list <- "chr1"
    window_width_list <- 1e6
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
  score_type_list <- c("corr_pearson", "dist_l1")
  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      chrom = chrom_list,
      window_width = window_width_list
    ),
    function(sample_type, chrom, window_width) {
      if (window_width >= 1e6) {
        for (strand in const.STRANDS) {
          line_graph.plot_line(
            sample_type = sample_type,
            chrom = chrom,
            strand = strand,
            ribo_nuc = const.RIBO_NUC_ALL,
            window_width = window_width,
            normalize_method = "divideNonGap",
            base_font_size = settings.BASE_FONT_SIZE
          )
        }
      }
      for (score_type in score_type_list) {
        for (strand in const.STRANDS) {
          line_graph.plot_matrix(
            sample_type = sample_type,
            chrom = chrom,
            strand = strand,
            ribo_nuc = const.RIBO_NUC_ALL,
            window_width = window_width,
            normalize_method = "divideNonGap",
            score_type = score_type,
            base_font_size = settings.BASE_FONT_SIZE
          )
        }
      }
      sample_type_enzyme_f <- paste0(sample_type, "_enzyme_f")
      if (window_width >= 1e5) {
        outlier_type_args <- list(
          p_value = list(list(p_value = 0.01), list(p_value = 0.05), list(p_value = 0.1)),
          sd_thresh = list(list(sd_thresh = 3), list(sd_thresh = 5), list(sd_thresh = 10))
        )
        for (outlier_type in names(outlier_type_args)) {
          for (outlier_args in outlier_type_args[[outlier_type]]) {
            for (strand in const.STRANDS) {
              line_graph.plot_outlier(
                sample_type = sample_type_enzyme_f,
                chrom = chrom,
                strand = strand,
                ribo_nuc = const.RIBO_NUC_ALL,
                window_width = window_width,
                normalize_method = "divideNonGap",
                outlier_type = outlier_type,
                outlier_args = outlier_args,
                base_font_size = settings.BASE_FONT_SIZE
              )
              line_graph.plot_outlier_indicators(
                sample_type = sample_type_enzyme_f,
                chrom = chrom,
                strand = strand,
                ribo_nuc = const.RIBO_NUC_ALL,
                window_width = window_width,
                normalize_method = "divideNonGap",
                outlier_type = outlier_type,
                outlier_args = outlier_args,
                base_font_size = settings.BASE_FONT_SIZE
              )
            }
          }
        }
      }
    }
  )
}
