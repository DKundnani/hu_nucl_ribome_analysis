dnaseq_analysis.plot_line <- function(
  sample_list,
  chrom,
  strand,
  window_width,
  assembly = "hg38",
  normalize_method = "divideNonGap",
  remove_outliers = TRUE,
  base_font_size = 11
) {
  if (chrom == const.CHROM_ALL) {
    stop("chrom cannot be", const.CHROM_ALL)
  }
  if (strand == const.STRAND_BOTH) {
    stop("strand cannot be", const.STRAND_BOTH)
  }
  data <- (
    purrr::map(
      sample_list,
      function(sample) {
        utils.load_dnaseq_binned(
          sample = sample,
          chrom = chrom,
          strand = strand,
          window_width = window_width,
          normalize_method = normalize_method,
          remove_outliers = remove_outliers,
          assembly = assembly
        ) |>
        tibble::as_tibble() |>
        dplyr::mutate(sample = !!sample, .before = 1) |>
        dplyr::rename(chrom = seqnames)
      }
    ) |>
    purrr::list_rbind()
  )

  tick_info <- utils.get_tick_info(utils.get_chrom_sizes(chrom))
  (
    ggplot2::ggplot(
      data = data,
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
    ggplot2::ylab("DNA-seq depth") +
    ggplot2::xlab("Chromosome position") +
    ggplot2::labs(
      title = stringr::str_c(
        "Binned DNA-seq coverage line graph\n",
        "Chromosome: ", chrom, "\n",
        "Assembly: ", assembly, "\n",
        "x-axis shows the position on the chromosome\n",
        "y-axis shows the coverage in the window\n",
        "Window width: ", utils.get_pretty_base_pairs(window_width), "\n",
        "Normalize method: ", normalize_method, "\n",
        "Outliers removed: ", remove_outliers
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
    ) +
    ggplot2::facet_wrap(
      dplyr::vars(sample),
      ncol = 1,
      scales = "free_y"
    )
  ) |>
  utils.write_ggplot(
    fn.dnaseq_analysis_plot(
      prefix = "line_graph",
      assembly = as.character(assembly),
      chrom = as.character(chrom),
      strand = as.character(strand),
      window_width = window_width,
      normalize_method = normalize_method,
      remove_outliers = remove_outliers
    ),
    width = 12,
    height = max(2 * length(sample_list), 2.5),
    limitsize = FALSE
  )
}

# Summary of the binned sequencing depth
dnaseq_analysis.get_summary_depth <- function(
  sample_list,
  chrom,
  strand,
  window_width,
  assembly = "hg38",
  normalize_method = "divideNonGap",
  remove_outliers = TRUE
) {
  data <- (
    purrr::map(
      sample_list,
      function(sample) {
        utils.load_dnaseq_binned(
          sample = sample,
          assembly = assembly,
          chrom = chrom,
          strand = strand,
          window_width = window_width,
          normalize_method = normalize_method,
          remove_outliers = remove_outliers
        ) |>
        utils.granges_to_tibble() |>
        dplyr::mutate(sample = !!sample, .before = 1)
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::mutate(sample = factor(sample, sample_list))
  )

  data_summary <- (
    data |>
    tidyr::drop_na() |>
    dplyr::group_by(sample) |>
    dplyr::summarize(
      n = dplyr::n(),
      depth_median = median(count),
      depth_mad = mad(count)
    )
  )

  # write tex version
  data_summary |>
  dplyr::mutate(
    depth_median = sprintf("%.2f", depth_median),
    depth_mad = sprintf("%.3f", depth_mad)
  ) |>
  xtable::xtable() |>
  utils.write_xtable(
    fn.dnaseq_analysis_data(
      ext = "tex",
      prefix = "summary_depth",
      assembly = as.character(assembly),
      chrom = as.character(chrom),
      strand = as.character(strand),
      window_width = window_width,
      normalize_method = normalize_method,
      remove_outliers = remove_outliers
    )
  )

  # write csv version
  data_summary |>
  utils.write_csv(
    fn.dnaseq_analysis_data(
      ext = "csv",
      prefix = "summary_depth",
      assembly = as.character(assembly),
      chrom = as.character(chrom),
      strand = as.character(strand),
      window_width = window_width,
      normalize_method = normalize_method,
      remove_outliers = remove_outliers
    )
  )

  # write outliers
  if (!remove_outliers) {
    data_outlier <- (
      data |>
      dplyr::filter(outlier) |>
      dplyr::mutate(
        start = as.integer(start),
        end = as.integer(end)
      )
    )

    data_outlier |>
    utils.write_csv(
      fn.dnaseq_analysis_data(
        ext = "csv",
        prefix = "summary_outlier",
        assembly = as.character(assembly),
        chrom = as.character(chrom),
        strand = as.character(strand),
        window_width = window_width,
        normalize_method = normalize_method,
        remove_outliers = remove_outliers
      )
    )

    data_outlier |>
    xtable::xtable() |>
    utils.write_xtable(
      fn.dnaseq_analysis_data(
        ext = "tex",
        prefix = "summary_outlier",
        assembly = as.character(assembly),
        chrom = as.character(chrom),
        strand = as.character(strand),
        window_width = window_width,
        normalize_method = normalize_method,
        remove_outliers = remove_outliers
      )
    )
  }
}

# Summary of the binary coverage (not binned)
dnaseq_analysis.get_summary_coverage <- function(
  sample_list,
  chrom,
  strand,
  assembly = "hg38",
  normalize_method = "divideNonGap",
  remove_outliers = TRUE
) {
  data <- (
    purrr::map(
      sample_list,
      function(sample) {
        utils.load_dnaseq_rds(
          sample = sample,
          chrom = chrom,
          assembly = assembly
        ) |>
        utils.granges_to_tibble() |>
        dplyr::filter(count > 0) |>
        dplyr::group_by(chrom, strand) |>
        dplyr::summarize(
          covered_size = sum(end - start + 1),
          .groups = "drop"
        ) |>
        dplyr::mutate(sample = !!sample, .before = 1)
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      chrom_size = utils.get_chrom_sizes(chrom),
      covered_frac = covered_size / chrom_size,
      uncovered_frac = 1 - covered_frac
    ) |>
    dplyr::select(
      sample,
      chrom,
      strand,
      covered_frac,
      uncovered_frac
    )
  )

  data |>
  utils.write_csv(
    fn.dnaseq_analysis_data(
      ext = "csv",
      prefix = "summary_coverage",
      assembly = as.character(assembly),
      chrom = as.character(chrom),
      normalize_method = normalize_method,
      remove_outliers = remove_outliers
    )
  )

  data |>
  dplyr::select(!c(chrom, strand)) |>
  xtable::xtable() |>
  utils.write_xtable(
    fn.dnaseq_analysis_data(
      ext = "tex",
      prefix = "summary_coverage",
      assembly = as.character(assembly),
      chrom = as.character(chrom),
      normalize_method = normalize_method,
      remove_outliers = remove_outliers
    )
  )
}

dnaseq_analysis.get_ribo_correlation_data <- function(
  dnaseq_list,
  ribo_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  cor_method = "spearman",
  normalize_method = "divideNonGap",
  remove_outliers = FALSE
) {
  if (length(intersect(dnaseq_list, ribo_list)) > 0) {
    stop("Overlapping DNA-seq/ribo names.")
  }

  data_dnaseq <- (
    purrr::map(
      dnaseq_list,
      function(sample) {
        utils.load_dnaseq_binned(
          sample = sample,
          chrom = chrom,
          strand = strand,
          window_width = window_width,
          assembly = "hg38",
          normalize_method = normalize_method,
          remove_outliers = remove_outliers
        ) |>
        utils.granges_to_tibble() |>
        dplyr::mutate(sample_dnaseq = !!sample, .before = 1)
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::select(
      sample_dnaseq,
      strand,
      start,
      end,
      count
    ) |>
    dplyr::nest_by(sample_dnaseq)
  )

  data_ribo <- (
    purrr::map(
      ribo_list,
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
        dplyr::mutate(sample_ribo = !!sample, .before = 1)
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::select(
      sample_ribo,
      strand,
      start,
      end,
      count
    ) |>
    dplyr::nest_by(sample_ribo)
  )

  dplyr::cross_join(
    data_dnaseq,
    data_ribo,
    suffix = c("_dnaseq", "_ribo")
  ) |>
  purrr::pmap(
    function(sample_dnaseq, sample_ribo, data_dnaseq, data_ribo) {
      utils.log(sample_dnaseq, sample_ribo)
      if (!all(data_dnaseq[["start"]] == data_ribo[["start"]])) {
        stop("Start positions do not match.")
      }
      tibble::tibble(
        chrom = !!chrom,
        strand = !!strand,
        sample_dnaseq = !!sample_dnaseq,
        sample_ribo = !!sample_ribo,
        cor = cor(
          data_dnaseq[["count"]],
          data_ribo[["count"]],
          method = cor_method,
          use = "complete.obs"
        )
      )
    }
  ) |>
  purrr::list_rbind()
}

dnaseq_analysis.plot_ribo_correlation <- function(
  dnaseq_list,
  sample_type_ribo,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  cor_method = "spearman",
  normalize_method = "divideNonGap",
  remove_outliers = FALSE,
  base_font_size = 11
) {
  cor_method_label <- c(
    spearman = "Spearman",
    pearson = "Pearson",
    kendall = "Kendall"
  )[[cor_method]]

  ribo_list <- const.get_sample_list(sample_type_ribo)
  data <- dnaseq_analysis.get_ribo_correlation_data(
    dnaseq_list = dnaseq_list,
    ribo_list = ribo_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width,
    cor_method = cor_method,
    normalize_method = normalize_method,
    remove_outliers = remove_outliers
  )

  (
    data |>
    dplyr::mutate(
      cor_label = sprintf("%.2f", cor),
      sample_dnaseq = const.get_sample_label_factor(sample_dnaseq, dnaseq_list),
      sample_ribo = const.get_sample_label_factor(sample_ribo, ribo_list)
    ) |>
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = sample_ribo,
        y = sample_dnaseq,
        fill = cor,
        label = cor_label
      )
    ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text() +
    ggplot2::labs(
      title = paste0(
        "DNA-seq coverage correlation matrix\n",
        "Ribo sample type = ", sample_type_ribo, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nucs = ", ribo_nuc, "\n",
        "Normalize method = ", normalize_method, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "X-axis shows the ribo sample\n",
        "Y-axis shows the DNA-seq sample\n",
        "Color shows the ", cor_method_label, " correlation"
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
      aspect.ratio = length(dnaseq_list) / length(ribo_list)
    ) +
    ggplot2::scale_fill_gradientn(
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      colors = c("blue", "white", "red"),
      na.value = "#fff0f0",
      oob = scales::squish
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = paste0(cor_method_label, " correlation"),
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
    fn.dnaseq_analysis_plot(
      sample_type_ribo = sample_type_ribo,
      prefix = "matrix",
      assembly = "hg38",
      chrom = chrom,
      strand = strand,
      window_width = window_width,
      ribo_nuc,
      cor_method,
      normalize_method = normalize_method,
      remove_outliers = remove_outliers
    ),
    width = max(length(ribo_list), 2) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = max(length(dnaseq_list), 4) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    limitsize = FALSE
  )
}

dnaseq_analysis.plot_line_ribo <- function(
  dnaseq_list,
  sample_type_ribo,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method_1 = "divideNonGap",
  normalize_method_2 = "mean",
  remove_outliers = TRUE,
  assembly = "hg38",
  base_font_size = 11
) {
  if (chrom == const.CHROM_ALL) {
    stop("Cannot use chromosome ", chrom, " for line plot")
  }
  if (strand == const.STRAND_BOTH) {
    stop("Cannot use strand ", strand, " for line plot")
  }

  data_dnaseq <- (
    purrr::map(
      dnaseq_list,
      function(sample) {
        utils.load_dnaseq_binned(
          sample = sample,
          chrom = chrom,
          strand = strand,
          window_width = window_width,
          assembly = "hg38",
          normalize_method = normalize_method_1,
          remove_outliers = remove_outliers
        ) |>
        utils.granges_to_tibble() |>
        dplyr::mutate(sample = !!sample, .before = 1) |>
        dplyr::mutate(count = utils.normalize(count, normalize_method_2)) |>
        dplyr::nest_by(sample, .keep = TRUE)
      }
    ) |>
    purrr::list_rbind()
  )

  ribo_list <- const.get_sample_list(sample_type_ribo)
  data_ribo <- (
    purrr::map(
      ribo_list,
      function(sample) {
        utils.load_ribos_binned(
          sample = sample,
          chrom = chrom,
          strand = strand,
          ribo_nuc = ribo_nuc,
          window_width = window_width,
          normalize_method = normalize_method_1
        ) |>
        utils.granges_to_tibble() |>
        dplyr::mutate(sample = !!sample, .before = 1) |>
        dplyr::mutate(count = utils.normalize(count, normalize_method_2)) |>
        dplyr::nest_by(sample, .keep = TRUE)
      }
    ) |>
    purrr::list_rbind()
  )

  geom_list <- (
    dplyr::cross_join(
      data_dnaseq,
      data_ribo,
      suffix = c("_dnaseq", "_ribo")
    ) |>
    purrr::pmap(
      function(sample_dnaseq, sample_ribo, data_dnaseq, data_ribo) {
        data_dnaseq <- (
          data_dnaseq |>
          dplyr::mutate(
            dnaseq = const.get_sample_label_factor(sample_dnaseq, dnaseq_list),
            ribo = const.get_sample_label_factor(sample_ribo, ribo_list),
            type = forcats::fct("dnaseq", c("dnaseq", "ribo")),
            row_key = forcats::fct_cross(type, ribo),
            row_col_key = forcats::fct_cross(row_key, dnaseq)
          )
        )
        data_ribo <- (
          data_ribo |>
          dplyr::mutate(
            dnaseq = const.get_sample_label_factor(sample_dnaseq, dnaseq_list),
            ribo = const.get_sample_label_factor(sample_ribo, ribo_list),
            type = forcats::fct("ribo", c("dnaseq", "ribo")),
            row_key = forcats::fct_cross(type, ribo),
            row_col_key = forcats::fct_cross(row_key, dnaseq)
          )
        )
        panel_dnaseq <- ggplot2::geom_line(
          data = data_dnaseq,
          mapping = ggplot2::aes(
            x = start,
            y = count,
            group = row_col_key,
            color = type
          ),
          linewidth = 0.5
        )
        panel_ribo <- ggplot2::geom_line(
          data = data_ribo,
          mapping = ggplot2::aes(
            x = start,
            y = count,
            group = row_col_key,
            color = type
          ),
          linewidth = 0.5
        )
        list(panel_dnaseq, panel_ribo)
      }
    ) |>
    purrr::list_c()
  )

  tick_info <- utils.get_tick_info(utils.get_chrom_sizes(chrom))
  (
    purrr::reduce(geom_list, `+`, .init = ggplot2::ggplot()) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_x_continuous(
      breaks = tick_info[["breaks"]],
      labels = tick_info[["labels"]]
    ) +
    ggplot2::ylab(
      if (normalize_method_2 == "none") {
        "Ribo frequency or sequencing depth"
      } else if (normalize_method_2 == "sum") {
        "Ribo frequency or sequencing depth (normalized by dividing by sum)"
      } else if (normalize_method_2 == "mean") {
        "Ribo frequency or sequencing depth (normalized by dividing by mean)"
      } else if (normalize_method_2 == "rank") {
        "Ribo frequency or sequencing depth (normalized by assigning rank)"
      } else {
        stop("Unknown normalize method: ", normalize_method_2)
      }
    ) +
    ggplot2::xlab("Chromosome position") +
    ggplot2::labs(
      title = paste0(
        "Ribo vs. DNA-seq\n",
        "Ribo sample type: ", sample_type_ribo, "\n",
        "Chromosome: ", chrom, "\n",
        "Strand: ", strand, "\n",
        "Ribo nucs: ", ribo_nuc, "\n",
        "X-axis shows the position on the chromosome\n",
        "Y-axis shows the normalized value of the window\n",
        "Normalize method 1: ", normalize_method_1, "\n",
        "Normalize method 2: ", normalize_method_2, "\n",
        "Window width: ", utils.get_pretty_base_pairs(window_width), "\n",
        "Outliers removed: ", remove_outliers
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
    ) +
    ggplot2::facet_grid(
      rows = dplyr::vars(row_key),
      cols = dplyr::vars(dnaseq),
      labeller = function(x) {
        if ("row_key" %in% colnames(x)) {
          x |>
          dplyr::rowwise() |>
          dplyr::mutate(
            row_key = if (stringr::str_starts(row_key, "dnaseq")) {
              "DNA-seq depth"
            } else if (stringr::str_starts(row_key, "ribo")) {
              row_key |>
              stringr::str_split(":") |>
              purrr::pluck(1, 2) # Get only the sample name
            } else {
              stop(stringr::str_c("Bad pattern: ", row_key))
            }
          ) |>
          dplyr::ungroup()
        } else {
          x
        }
      },
      scales = "free_y"
    )
  ) |>
  utils.write_ggplot(
    fn.dnaseq_analysis_plot(
      sample_type_ribo = sample_type_ribo,
      prefix = "line_graph_ribo",
      assembly = "hg38",
      chrom = as.character(chrom),
      strand = as.character(strand),
      window_width = window_width,
      ribo_nuc,
      normalize_method_2,
      normalize_method = normalize_method_1,
      remove_outliers = remove_outliers
    ),
    width = 12 * length(dnaseq_list) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = max(4 * length(ribo_list), 6),
    limitsize = FALSE
  )
}

# Check how ribonucleotides fall in covered or uncovered regions
dnaseq_analysis.get_ribo_coverage <- function(
  dnaseq_list,
  sample_type_ribo,
  chrom,
  assembly = "hg38"
) {
  data_dnaseq <- (
    purrr::map(
      dnaseq_list,
      function(sample) {
        ranges <- utils.load_dnaseq_rds(
          sample = sample,
          chrom = chrom
        )
        ranges <- ranges[ranges$count > 0]
        ranges <- GenomicRanges::reduce(ranges)
        tibble::tibble(
          sample = !!sample,
          ranges = list(!!ranges)
        )
      }
    ) |>
    purrr::list_rbind()
  )

  ribo_list <- const.get_sample_list(sample_type_ribo)
  data_ribo <- (
    purrr::map(
      ribo_list,
      function(sample) {
        ranges <- utils.load_ribos_rds(
          sample = sample,
          chrom = chrom
        )
        tibble::tibble(
          sample = !!sample,
          ranges = list(!!ranges)
        )
      }
    ) |>
    purrr::list_rbind()
  )

  table_coverage <- (
    dplyr::cross_join(
      data_dnaseq,
      data_ribo,
      suffix = c("_dnaseq", "_ribo")
    ) |>
    purrr::pmap(
      function(
        sample_dnaseq,
        sample_ribo,
        ranges_dnaseq,
        ranges_ribo
      ) {
        ribos = NROW(ranges_ribo)
        ribos_covered <- sum(ranges_ribo %over% ranges_dnaseq)
        ribos_uncovered <- ribos - ribos_covered
        ribos_uncovered_frac <- ribos_uncovered / ribos
        tibble::tibble(
          chrom,
          sample_dnaseq,
          sample_ribo,
          ribos,
          ribos_covered,
          ribos_uncovered,
          ribos_uncovered_frac
        )
      }
    ) |>
    purrr::list_rbind()
  )

  table_coverage |>
  dplyr::mutate(
    sample_dnaseq = const.get_sample_label_factor(sample_dnaseq, dnaseq_list),
    sample_ribo = const.get_sample_label_factor(sample_ribo, ribo_list)
  ) |>
  dplyr::arrange(sample_dnaseq, sample_ribo) |>
  xtable::xtable() |>
  utils.write_xtable(
    fn.dnaseq_analysis_data(
      ext = "tex",
      sample_type_ribo = sample_type_ribo,
      prefix = "ribo_coverage",
      assembly = assembly,
      chrom = chrom,
      remove_outliers = FALSE
    )
  )

  table_coverage |>
  utils.write_csv(
    fn.dnaseq_analysis_data(
      ext = "csv",
      sample_type_ribo = sample_type_ribo,
      prefix = "ribo_coverage",
      assembly = assembly,
      chrom = chrom,
      sample_type_ribo,
      remove_outliers = FALSE
    )
  )
}

dnaseq_analysis.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_ribo_list <- c(const.SAMPLE_TEST)
    chrom_list <- "chr1"
    window_width_list <- c(1e5, 1e6)
  } else {
    sample_type_ribo_list <- c(const.SAMPLE_NORMAL)
    chrom_list <- const.CHROMS
    window_width_list <- c(1e3, 1e4, 1e5, 1e6)
  }

  for (chrom in chrom_list) {
    if (test) {
      dnaseq_analysis.get_ribo_coverage(
        dnaseq_list = "F",
        sample_type_ribo = const.SAMPLE_TEST_ENZYME_F,
        chrom = chrom
      )
    } else {
      dnaseq_analysis.get_ribo_coverage(
        dnaseq_list = "F",
        sample_type_ribo = const.SAMPLE_NORMAL_ENZYME_F,
        chrom = chrom
      )
    }
    dnaseq_analysis.get_summary_coverage(
      sample_list = const.SAMPLES_DNASEQ,
      chrom = chrom,
      remove_outliers = FALSE
    )
    for (window_width in window_width_list) {
      dnaseq_analysis.get_summary_depth(
        sample_list = const.SAMPLES_DNASEQ,
        chrom = "chr1",
        strand = const.STRAND_BOTH,
        window_width = window_width,
        remove_outliers = FALSE
      )

      for (sample_type_ribo in sample_type_ribo_list) {
        dnaseq_analysis.plot_ribo_correlation(
          dnaseq_list = "F",
          sample_type_ribo = paste0(sample_type_ribo, "_enzyme_f"),
          chrom = const.CHROM_ALL,
          strand = const.STRAND_BOTH,
          ribo_nuc = const.RIBO_NUC_ALL,
          window_width = window_width,
          cor_method = "spearman",
          remove_outliers = FALSE,
          base_font_size = settings.BASE_FONT_SIZE
        )
      }

      if (window_width >= 1e6) {
        for (strand in const.STRANDS) {
          dnaseq_analysis.plot_line(
            sample_list = const.SAMPLES_DNASEQ,
            window_width = window_width,
            chrom = chrom,
            strand = strand,
            remove_outliers = TRUE,
            base_font_size = settings.BASE_FONT_SIZE
          )
          for (sample_type_ribo in sample_type_ribo_list) {
            dnaseq_analysis.plot_line_ribo(
              dnaseq_list = "F",
              sample_type_ribo = paste0(sample_type_ribo, "_enzyme_f"),
              chrom = chrom,
              strand = strand,
              ribo_nuc = const.RIBO_NUC_ALL,
              window_width = window_width,
              normalize_method_1 = "divideNonGap",
              normalize_method_2 = "mean",
              remove_outliers = TRUE,
              base_font_size = settings.BASE_FONT_SIZE
            )
          }
        }
      }
    }
  }
}
