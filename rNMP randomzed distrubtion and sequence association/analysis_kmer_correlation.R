kmer_correlation.load_kmer_data <- function(
  kmer_size,
  kmer_list,
  chrom,
  strand,
  window_width
) {
  utils.log(
    "Loading kmer data:",
    chrom,
    strand,
    utils.get_pretty_base_pairs(window_width),
    kmer_size,
    paste0("(", paste0(kmer_list, collapse = ", "), ")")
  )

  chrom_list <- const.expand_chrom(chrom)
  strand_list <- const.expand_strand(strand)

  purrr::map(
    chrom_list,
    function(chrom) {
      ranges <- (
        fn.kmer_binned(
          chrom = chrom,
          window_width = window_width,
          kmer_size = kmer_size
        ) |>
        utils.read_rds() |>
        utils.ignore_log() # Do not log the file name
      )
      ranges <- ranges[GenomicRanges::strand(ranges) %in% strand_list, ]
      count_total <- rowSums(ranges$count)
      ranges$count <- sweep(ranges$count, 1, count_total, `/`) # convert to frequencies
      ranges$count <- ranges$count[, kmer_list]
      
      dplyr::bind_cols(
        tibble::tibble(
          chrom = !!chrom,
          strand = as.vector(GenomicRanges::strand(ranges)),
          start = as.vector(GenomicRanges::start(ranges)),
          end = as.vector(GenomicRanges::end(ranges))
        ),
        ranges$count
      ) |>
      tidyr::pivot_longer(
        dplyr::all_of(kmer_list),
        names_to = "kmer",
        values_to = "count"
      ) |>
      dplyr::select(
        kmer,
        chrom,
        strand,
        start,
        end,
        count
      )
    }
  ) |>
  purrr::list_rbind() |>
  dplyr::arrange(kmer, chrom, strand, start)
}

kmer_correlation.load_ribo_data <- function(
  sample_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method = "divideNonGap"
) {
  chrom_list <- const.expand_chrom(chrom)
  strand_list <- const.expand_strand(strand)
  ribo_nuc_list <- const.expand_ribo_nuc(ribo_nuc)

  purrr::map(
    chrom_list,
    function(chrom) {
      weights <- (
        fn.weights_binned(chrom = chrom, window_width = window_width) |>
        utils.read_rds()
      )

      purrr::map(
        sample_list,
        function(sample) {
          ranges <- (
            fn.ribos_binned(
              sample = sample,
              chrom = chrom,
              window_width = window_width
            ) |>
            utils.read_rds()
          )
          if (is.null(ranges)) {
            # Means there was a missing file
            ranges <- utils.get_ranges_binned(
              chrom = chrom,
              window_width = window_width,
              strands = TRUE
            )
            ranges$count <- matrix(0, nrow = NROW(ranges), ncol = length(ribo_nuc_list))
            colnames(ranges$count) <- ribo_nuc_list
          }
          ranges <- ranges[GenomicRanges::strand(ranges) %in% strand_list]
          ranges$count <- ranges$count[, ribo_nuc_list] |> rowSums()
          if (!all(GenomicRanges::start(ranges) == GenomicRanges::start(weights))) {
            stop("Start positions of ribos and weights do not match")
          }
          if (normalize_method == "divideNonGap") {
            ranges$count <- ranges$count / weights$count
          } else if (normalize_method == "none") {
          } else {
            stop("Invalid normalize_method")
          }
          ranges$count[!is.finite(ranges$count)] <- NA
          tibble::tibble(
            sample = !!sample,
            chrom = !!chrom,
            strand = as.vector(GenomicRanges::strand(ranges)),
            start = as.vector(GenomicRanges::start(ranges)),
            end = as.vector(GenomicRanges::end(ranges)),
            count = as.vector(ranges$count)
          )
        }
      ) |>
      purrr::list_rbind()
    }
  ) |>
  purrr::list_rbind() |>
  dplyr::arrange(sample, chrom, strand, start)
}

kmer_correlation.load_positional_data <- function(
  sample_list,
  kmer_size,
  kmer_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  separate_ribo_nucs = FALSE,
  normalize_method_1 = "divideNonGap",
  normalize_method_2 = "mean"
) {
  data_kmer <- (
    kmer_correlation.load_kmer_data(
      kmer_size = kmer_size,
      kmer_list = kmer_list,
      chrom = chrom,
      strand = strand,
      window_width = window_width
    ) |>
    dplyr::group_by(kmer) |>
    dplyr::mutate(count = utils.normalize(count, normalize_method_2)) |>
    dplyr::ungroup()
  )
  data_ribo <- (
    kmer_correlation.load_ribo_data(
      sample_list = sample_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_width = window_width,
      normalize_method = normalize_method_1
    ) |>
    dplyr::group_by(sample) |>
    dplyr::mutate(count = utils.normalize(count, normalize_method_2)) |>
    dplyr::ungroup()
  )
  list(kmer = data_kmer, ribo = data_ribo)
}

kmer_correlation.load_correlation_data <- function(
  sample_list,
  kmer_size,
  kmer_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  cor_method = "spearman",
  normalize_method = "divideNonGap"
) {
  data <- kmer_correlation.load_positional_data(
    sample_list = sample_list,
    kmer_size = kmer_size,
    kmer_list = kmer_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width
  )

  chrom_list <- const.expand_chrom(chrom)
  strand_list <- const.expand_strand(strand)

  data[["ribo"]] <- (
    data[["ribo"]] |>
    dplyr::mutate(
      chrom = factor(chrom, chrom_list),
      strand = factor(strand, strand_list)
    ) |>
    dplyr::nest_by(sample) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data = (
        data |>
        dplyr::arrange(chrom, strand, start) |>
        list()
      )
    ) |>
    dplyr::ungroup()
  )

  data[["kmer"]] <- (
    data[["kmer"]] |>
    dplyr::mutate(
      chrom = factor(chrom, chrom_list),
      strand = factor(strand, strand_list)
    ) |>
    dplyr::nest_by(kmer) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      data = (
        data |>
        dplyr::arrange(chrom, strand, start) |>
        list()
      )
    ) |>
    dplyr::ungroup()
  )

  dplyr::cross_join(
    data[["ribo"]],
    data[["kmer"]],
    suffix = c("_ribo", "_kmer")
  ) |>
  purrr::pmap(
    function(sample, kmer, data_ribo, data_kmer) {
      utils.log("Correlation:", sample, kmer)

      if (!all(data_ribo[["chrom"]] == data_kmer[["chrom"]])) {
        stop("Chromosomes don't match")
      }
      if (!all(data_ribo[["strand"]] == data_kmer[["strand"]])) {
        stop("Strands don't match")
      }
      if (!all(data_ribo[["start"]] == data_kmer[["start"]])) {
        stop("Starts don't match")
      }
      if (!all(is.na(data_ribo[["count"]]) == is.na(data_kmer[["count"]]))) {
        stop("NAs don't match")
      }
      
      tibble::tibble(
        sample = !!sample,
        kmer = !!kmer,
        cor = cor(
          data_ribo[["count"]],
          data_kmer[["count"]],
          method = cor_method,
          use = "complete.obs"
        )
      )
    }
  ) |>
  purrr::list_rbind()
}

kmer_correlation.plot_correlation <- function(
  sample_type,
  kmer_size,
  kmer_list = NULL,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  cor_method = "spearman",
  normalize_method = "divideNonGap",
  base_font_size = 11,
  cor_font_size = 8, # Font size of the correlation text.
  rev_kmer_list = FALSE,
  aspect_ratio = NULL, # Aspect ratio of the grid; if NULL, it is automatically determined.
  remove_axis_labels = FALSE, # Whether to remove the axis labels.
  # Colors must be in the order of -1, 0, 1 correlation.
  # Default colors: purple, gray, green.
  colors = c("#af8dc3", "#f7f7f7", "#7fbf7b"),
  use_pheatmap = TRUE # Whether to use pheatmap instead of ggplot2.
) {
  sample_list <- const.get_sample_list(sample_type)
  if (is.null(kmer_list)) {
    kmer_list <- kmer.get_all_kmers(kmer_size)
  }
  if (is.null(aspect_ratio)) {
    aspect_ratio <- length(kmer_list) / length(sample_list)
  }
  data <- kmer_correlation.load_correlation_data(
    sample_list = sample_list,
    kmer_size = kmer_size,
    kmer_list = kmer_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width,
    cor_method = cor_method,
    normalize_method = normalize_method
  ) |>
  dplyr::mutate(
    cor_label = sprintf("%.2f", cor),
    kmer = factor(kmer, if (rev_kmer_list) rev(kmer_list) else kmer_list)
  )

  file_out <- fn.kmer_correlation_plot(
    sample_type = sample_type,
    prefix = "matrix",
    window_width = window_width,
    kmer_size = kmer_size,
    if (length(kmer_list) == kmer.get_num_kmers(kmer_size)) {
      "all"
    } else {
      paste0(kmer_list, collapse = "_")
    },
    chrom,
    strand,
    ribo_nuc,
    cor_method,
    normalize_method
  )

  if (use_pheatmap) {
    library(pheatmap)
    data <- (
      data |>
      tidyr::pivot_wider(
        names_from = sample,
        values_from = cor,
        id_cols = kmer
      ) |>
      dplyr::select(kmer, dplyr::all_of(sample_list)) |>
      dplyr::arrange(dplyr::desc(kmer))
    )
    mat <- as.matrix(data |> dplyr::select(!kmer))
    rownames(mat) <- data[["kmer"]]
    colnames(mat) <- (
      colnames(data |> dplyr::select(!kmer)) |>
      purrr::map_chr(const.get_orig_sample)
    )

    utils.log("(output)", file_out)
    utils.create_parent_dir(file_out)
    png(
      file_out,
      width = 1 + (ncol(mat) * 0.4),
      height = 1 + (nrow(mat) * 0.35),
      unit = "in",
      res = 600
    )
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      breaks = seq(-1, 1, length.out = 200),
      fontsize_number = 12,
      number_color = "black",
      fontsize_row = 12,
      fontsize_col = 12,
      color = grDevices::colorRampPalette(colors)(200),
      labels_row = rownames(mat),
      labels_col = colnames(mat)
    )
    dev.off()
  } else {
    (
      data |>
      dplyr::mutate(
        sample = const.get_sample_label_factor(sample, sample_list)
      ) |>
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = sample,
          y = kmer,
          fill = cor,
          label = cor_label
        )
      ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(size = cor_font_size) +
      ggplot2::labs(
        title = paste0(
          "Kmer correlation matrix\n",
          "Sample type: ", sample_type, "\n",
          "Chromosome: ", chrom, "\n",
          "Strand: ", strand, "\n",
          "Ribo nuc: ", ribo_nuc, "\n",
          "Window width: ", utils.get_pretty_base_pairs(window_width), "\n",
          "Normalize method: ", normalize_method, "\n",
          "Correlation method: ", cor_method, "\n",
          "Kmer size: ", utils.get_pretty_base_pairs(kmer_size), "\n",
          "X-axis shows the kmer\n",
          "Y-axis shows the sample\n",
          "Color shows correlation\n"
        )
      ) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("Kmer") +
      ggplot2::theme_classic(base_size = base_font_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        legend.position = "right",
        axis.title = if (remove_axis_labels) {
          ggplot2::element_blank()
        } else {
          ggplot2::element_text()
        },
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        aspect.ratio = aspect_ratio
      ) +
      ggplot2::scale_fill_gradientn(
        name = paste0(
          stringr::str_to_sentence(cor_method),
          " correlation"
        ),
        limits = c(-1, 1),
        breaks = c(-1, 0, 1),
        colors = colors,
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
      file_out,
      width = max(length(sample_list), 2),
      height = (
        (max(length(sample_list), 2) * aspect_ratio) +
        settings.PLOT_SAMPLE_LABEL_MARGIN
      ),
      limitsize = FALSE
    )
  }
}

kmer_correlation.plot_line <- function(
  sample_type,
  kmer_size,
  kmer_list = NULL,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method_1 = "divideNonGap",
  normalize_method_2 = "mean",
  base_font_size = 11
) {
  if (chrom == const.CHROM_ALL) {
    stop("Cannot plot line graph of all chromosomes")
  }
  if (strand == const.STRAND_BOTH) {
    stop("Cannot plot line graph of both strands")
  }
  sample_list = const.get_sample_list(sample_type)
  if (is.null(kmer_list)) {
    kmer_list <- kmer.get_all_kmers(kmer_size)
  }
  data <- kmer_correlation.load_positional_data(
    sample_list = sample_list,
    kmer_size = kmer_size,
    kmer_list = kmer_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width,
    normalize_method_1 = normalize_method_1,
    normalize_method_2 = normalize_method_2
  )

  tick_info <- utils.get_tick_info(utils.get_chrom_sizes(chrom))

  (
    ggplot2::ggplot() +
    ggplot2::geom_line(
      data = (
        data[["kmer"]] |>
        tidyr::expand_grid(sample = sample_list) |>
        dplyr::mutate(
          sample = const.get_sample_label_factor(sample, sample_list),
          kmer = forcats::fct(kmer, kmer_list),
          type = forcats::fct("kmer", c("kmer", "ribo")),
          row_key = forcats::fct_cross(type, sample, keep_empty = TRUE),
          row_col_key = forcats::fct_cross(row_key, kmer, keep_empty = TRUE)
        )
      ),
      mapping = ggplot2::aes(
        x = start,
        y = count,
        group = row_col_key,
        color = type
      )
    ) +
    ggplot2::geom_line(
      data = (
        data[["ribo"]] |>
        tidyr::expand_grid(kmer = kmer_list) |>
        dplyr::mutate(
          sample = const.get_sample_label_factor(sample, sample_list),
          kmer = forcats::fct(kmer, kmer_list),
          type = forcats::fct("ribo", c("kmer", "ribo")),
          row_key = forcats::fct_cross(type, sample, keep_empty = TRUE),
          row_col_key = forcats::fct_cross(row_key, kmer, keep_empty = TRUE)
        )
      ),
      mapping = ggplot2::aes(
        x = start,
        y = count,
        group = row_col_key,
        color = type
      )
    ) +
    ggplot2::geom_line(linewidth = 0.5) +
    (
      if ((normalize_method_2 == "mean") && (sample_type == const.SAMPLE_DNASEQ)) {
        ggplot2::coord_cartesian(ylim = c(0, 2))
      } else {
        ggplot2::geom_blank()
      }
    ) +
    ggplot2::scale_x_continuous(
      breaks = tick_info[["breaks"]],
      labels = tick_info[["labels"]]
    ) +
    ggplot2::ylab(
      if (normalize_method_2 == "mean") {
        "Enrichment factor (count / [mean count])"
      } else if (normalize_method_2 == "sum") {
        "Density (count / [sum count])"
      } else if (normalize_method_2 == "rank") {
        "Count rank"
      } else if (normalize_method_2 == "none") {
        "Count"
      } else {
        stop("Invalid normalize_method_2: ", normalize_method_2)
      }
    ) +
    ggplot2::xlab("Chromosome position") +
    ggplot2::labs(
      title = paste0(
        "Ribo and kmer vs. chromsome position\n",
        "Sample type = ", sample_type, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nucs = ", ribo_nuc, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "Kmer size = ", utils.get_pretty_base_pairs(kmer_size), "\n",
        "Normalize method 1 = ", normalize_method_1, "\n",
        "Normalize method 2 = ", normalize_method_2, "\n",
        "X-axis shows the position on the chromosome\n",
        "Y-axis shows the value in the window"
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
      cols = dplyr::vars(kmer),
      labeller = function(x) {
        if ("row_key" %in% colnames(x)) {
          x |>
          dplyr::rowwise() |>
          dplyr::mutate(
            row_key = if (stringr::str_starts(row_key, "kmer")) {
              "kmer"
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
    fn.kmer_correlation_plot(
      sample_type = sample_type,
      prefix = "line",
      window_width = window_width,
      kmer_size = kmer_size,
      (
        if (length(kmer_list) == kmer.get_num_kmers(kmer_size)) {
          "all"
        } else {
          paste0(kmer_list, collapse = "_")
        }
      ),
      chrom,
      strand,
      ribo_nuc,
      normalize_method_1,
      normalize_method_2
    ),
    width = 12 * length(kmer_list) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = max(3 * length(sample_list), 6),
    limitsize = FALSE
  )
}

kmer_correlation.plot_scatter <- function(
  sample_type,
  kmer_size,
  kmer_list = NULL,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method_1 = "divideNonGap",
  normalize_method_2 = "mean",
  base_font_size = 11
) {
  sample_list = const.get_sample_list(sample_type)
  if (is.null(kmer_list)) {
    kmer_list <- kmer.get_all_kmers(kmer_size)
  }
  data <- kmer_correlation.load_positional_data(
    sample_list = sample_list,
    kmer_size = kmer_size,
    kmer_list = kmer_list,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_width = window_width,
    normalize_method_1 = normalize_method_1,
    normalize_method_2 = normalize_method_2
  )

  normalize_label <- if (normalize_method_2 == "mean") {
    "Enrichment factor (count / [count mean])"
  } else if (normalize_method_2 == "sum") {
    "Density (count / [count sum])"
  } else if (normalize_method_2 == "rank") {
    "Count rank"
  } else if (normalize_method_2 == "none") {
    "Count"
  } else {
    stop("Invalid normalize_method_2: ", normalize_method_2)
  }

  (
    dplyr::cross_join(
      (
        data[["kmer"]] |>
        dplyr::select(kmer, chrom, strand, start, count) |>
        dplyr::nest_by(kmer, .keep = TRUE, .key = "data_kmer")
      ),
      (
        data[["ribo"]] |>
        dplyr::select(sample, chrom, strand, start, count) |>
        dplyr::nest_by(sample, .keep = TRUE, .key = "data_ribo")
      )
    ) |>
    dplyr::rowwise() |>
    dplyr::transmute(
      data_both = (
        dplyr::inner_join(
          data_kmer,
          data_ribo,
          by = c("chrom", "strand", "start"),
          suffix = c("_kmer", "_ribo")
        ) |>
        list()
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::pull(data_both) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      sample = const.get_sample_label_factor(sample, sample_list),
      kmer = factor(kmer, kmer_list)
    ) |>
    ggplot2::ggplot(mapping = ggplot2::aes(count_kmer, count_ribo)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::xlab(paste0("Kmer ", normalize_label)) +
    ggplot2::ylab(paste0("Ribo ", normalize_label)) +
    ggplot2::labs(
      title = paste0(
        "Ribo vs. kmer scatter plot (each dot is one window)\n",
        "Sample type = ", sample_type, "\n",
        "Chromosome = ", chrom, "\n",
        "Strand = ", strand, "\n",
        "Ribo nuc = ", ribo_nuc, "\n",
        "Window width = ", utils.get_pretty_base_pairs(window_width), "\n",
        "Kmer size = ", utils.get_pretty_base_pairs(kmer_size), "\n",
        "Normalize method 1 = ", normalize_method_1, "\n",
        "Normalize method 2 = ", normalize_method_2, "\n",
        "X-axis shows the ribo values\n",
        "Y-axis shows the kmer values"
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
      rows = dplyr::vars(sample),
      cols = dplyr::vars(kmer),
      scales = "free"
    )
  ) |>
  utils.write_ggplot(
    fn.kmer_correlation_plot(
      sample_type = sample_type,
      prefix = "scatter",
      window_width = window_width,
      kmer_size = kmer_size,
      (
        if (length(kmer_list) == kmer.get_num_kmers(kmer_size)) {
          "all"
        } else {
          paste0(kmer_list, collapse = "_")
        }
      ),
      chrom,
      strand,
      ribo_nuc,
      normalize_method_1,
      normalize_method_2
    ),
    width = max(4 * length(kmer_list), 6) + settings.PLOT_SAMPLE_LABEL_MARGIN,
    height = max(4 * length(sample_list), 6),
    limitsize = FALSE
  )
}

kmer_correlation.do_main <- function(test = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    chrom_list_corr <- "chr1"
    chrom_list_line <- "chr1"
    window_width_list <- 1e6
    kmer_size_list <- c(1, 2)
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ_DRAW
    )
    chrom_list_corr <- const.CHROMS
    chrom_list_line <- const.CHROMS
    window_width_list <- c(1e3, 1e4, 1e5, 1e6)
    kmer_size_list <- c(1, 2)
  }

  tidyr::expand_grid(
    sample_type = sample_type_list,
    chrom = chrom_list_corr,
    strand = const.STRAND_BOTH,
    ribo_nuc = const.RIBO_NUC_ALL
  ) |>
  purrr::pwalk(
    function(sample_type, chrom, strand, ribo_nuc) {
      tidyr::expand_grid(
        window_width = window_width_list,
        kmer_size = kmer_size_list
      ) |>
      purrr::pwalk(
        function(window_width, kmer_size) {
          kmer_correlation.plot_correlation(
            sample_type = sample_type,
            kmer_size = kmer_size,
            chrom = chrom,
            strand = strand,
            ribo_nuc = ribo_nuc,
            window_width = window_width,
            cor_method = "pearson",
            normalize_method = "divideNonGap",
            base_font_size = settings.BASE_FONT_SIZE
          )
        }
      )
    }
  )

  tidyr::expand_grid(
    sample_type = sample_type_list,
    window_width = 1e6,
    chrom = chrom_list_line,
    strand = const.STRANDS,
    ribo_nuc = const.RIBO_NUC_ALL
  ) |>
  purrr::pwalk(
    function(sample_type, window_width, chrom, strand, ribo_nuc) {
      kmer_correlation.plot_line(
        sample_type = paste0(sample_type, "_enzyme_f"),
        kmer_size = 1,
        chrom = chrom,
        strand = strand,
        ribo_nuc = ribo_nuc,
        window_width = window_width,
        normalize_method_1 = "divideNonGap",
        normalize_method_2 = "mean",
        base_font_size = settings.BASE_FONT_SIZE
      )
      kmer_correlation.plot_scatter(
        sample_type = paste0(sample_type, "_enzyme_f"),
        kmer_size = 1,
        chrom = chrom,
        strand = strand,
        ribo_nuc = ribo_nuc,
        window_width = window_width,
        normalize_method_1 = "divideNonGap",
        normalize_method_2 = "mean",
        base_font_size = settings.BASE_FONT_SIZE
      )
    }
  )
}
