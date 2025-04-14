library(BSgenome)
library(plyranges)

kmer_heatmaps.save_kmer_data <- function(
  sample_list,
  window_radius,
  kmer_size,
  overwrite = FALSE
) {
  purrr::pwalk(
    tidyr::expand_grid(
      sample = sample_list,
      chrom = const.CHROMS
    ),
    function(sample, chrom) {
      utils.log(sample, chrom, window_radius, kmer_size)

      file_out <- fn.kmer_heatmap_data(
        sample = sample,
        chrom = chrom,
        window_radius = window_radius,
        kmer_size = kmer_size
      )
      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }

      granges <- utils.load_ribos_rds(sample, chrom)
      if (is.null(granges) || (NROW(granges) == 0)) {
        utils.log("No ribos for", sample, chrom)
        utils.write_rds(NULL, file_out)
        return(invisible())
      }
      purrr::map(
        c("up", "down"),
        function(window_side) {
          (
            if (window_side == "up") {
              plyranges::flank_upstream(granges, width = window_radius)
            } else if (window_side == "down") {
              plyranges::flank_downstream(granges, width = window_radius)
            } else {
              stop("Invalid window_side:", window_side)
            }
          ) |>
          granges => dplyr::bind_cols(
            strand = as.vector(GenomicRanges::strand(granges)),
            ribo_nuc = granges$ribo_nuc,
            count = granges$count,
            count_kmer = (
              BSgenome::BSgenomeViews(
                BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                granges
              ) |>
              Biostrings::oligonucleotideFrequency(width = kmer_size)
            )
          ) |>
          dplyr::nest_by(strand, ribo_nuc) |>
          purrr::pmap(
            function(strand, ribo_nuc, data) {
              count <- data$count
              count_kmer <- data$count_kmer

              kmer_list <- colnames(count_kmer)
              count_kmer <- sweep(count_kmer, 1, count, "*") # multiply each row by corresponding count
              count_kmer <- colSums(count_kmer) # sum across rows
              count_kmer <- tibble::as_tibble_row(count_kmer)
              colnames(count_kmer) <- kmer_list

              count <- sum(count) # sum across rows

              tibble::tibble(
                chrom = !!chrom,
                strand = !!strand,
                ribo_nuc = !!ribo_nuc,
                window_side = !!window_side,
                count = !!count,
                count_kmer = !!count_kmer
              )
            }
          ) |>
          purrr::list_rbind()
        }
      ) |>
      purrr::list_rbind() |>
      dplyr::arrange(chrom, strand, ribo_nuc, window_side) |>
      utils.write_rds(file_out)
    }
  )
}

kmer_heatmaps.save_kmer_data_bg <- function(kmer_size, overwrite = FALSE) {
  purrr::walk(
    const.CHROMS,
    function(chrom) {
      file_out <- fn.kmer_heatmap_data_bg(
        chrom = chrom,
        kmer_size = kmer_size
      )

      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }

      GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(
          start = 1,
          end = utils.get_chrom_sizes(chrom)
        ),
        strand = "+"
      ) |>
      granges => {
        count_plus <- (
          BSgenome::BSgenomeViews(
            BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
            granges
          ) |>
          Biostrings::oligonucleotideFrequency(width = kmer_size) |>
          tibble::as_tibble()
        )
        count_minus <- count_plus
        colnames(count_minus) <- kmer.reverse_complement(colnames(count_minus))
        count_minus <- count_minus[, colnames(count_plus)]
        
        dplyr::bind_rows(
          tibble::tibble(
            chrom = as.vector(GenomicRanges::seqnames(granges)),
            strand = "+",
            count_kmer = count_plus
          ),
          tibble::tibble(
            chrom = as.vector(GenomicRanges::seqnames(granges)),
            strand = "-",
            count_kmer = count_minus
          )
        )
      } |>
      dplyr::arrange(chrom, strand) |>
      utils.write_rds(file_out)
    }
  )
}

kmer_heatmaps.load_data <- function(
  sample_list,
  chrom,
  strand,
  ribo_nuc,
  window_radius,
  window_side,
  kmer_size,
  normalize_method = "ratio" # "none", "ratio", "sum1"
) {
  chrom_list <- const.expand_chrom(chrom)
  strand_list <- const.expand_strand(strand)
  ribo_nuc_list <- const.expand_ribo_nuc(ribo_nuc)
  window_side_list <- if (window_side == "both") {
    c("up", "down")
  } else if (window_side %in% c("up", "down")){
    window_side
  } else {
    stop(paste0("Unknown window side: ", window_side))
  }

  utils.log("Load BG")
  data_bg <- (
    purrr::map(
      chrom_list,
      function(chrom) {
        fn.kmer_heatmap_data_bg(
          chrom = chrom,
          kmer_size = kmer_size
        ) |>
        utils.read_rds() |>
        dplyr::filter(strand %in% strand_list) |>
        utils.ignore_log()
      }
    ) |>
    purrr::list_rbind() |>
    x => colSums(x$count_kmer) |> # sum over the windows
    x => (x / sum(x)) # convert to frequencies
  )

  purrr::map(
    sample_list,
    function(sample) {
      utils.log("Load", sample)

      data <- (
        purrr::map(
          chrom_list,
          function(chrom) {
            fn.kmer_heatmap_data(
              sample = sample,
              chrom = chrom,
              window_radius = window_radius,
              kmer_size = kmer_size
            ) |>
            utils.read_rds() |>
            x => (
              if (is.null(x)) {
                x
              } else {
                x |>
                dplyr::filter(
                  (strand %in% strand_list) &
                  (ribo_nuc %in% ribo_nuc_list) &
                  (window_side %in% window_side_list)
                )
              }
            ) |>
            utils.ignore_log()
          }
        ) |>
        purrr::list_rbind() |>
        x => colSums(x$count_kmer) |> # sum over the windows
        x => (x / sum(x)) # convert to frequencies
      )
    
      if (normalize_method == "ratio") {
        data <- data / data_bg
      } else if (normalize_method == "sum1") {
        data <- data / data_bg
        data <- data / sum(data)
      } else if (normalize_method == "none") {
      } else {
        stop(paste0("Unknown normalize method: ", normalize_method))
      }

      tibble::tibble(
        sample = !!sample,
        kmer = names(data),
        value = data
      )
    }
  ) |>
  purrr::list_rbind()
}

kmer_heatmaps.plot <- function(
  sample_type,
  chrom = const.CHROM_ALL,
  strand = const.STRAND_BOTH,
  ribo_nuc = const.RIBO_NUC_ALL,
  window_side = "both",
  window_radius = 50,
  kmer_size = 1,
  normalize_method = "ratio", # "none", "ratio", or "sum1"
  base_font_size = 11,
  freq_font_size = 8, # Font size of the frequency text.
  rev_kmer_list = FALSE,
  aspect_ratio = NULL, # Aspect ratio of the grid; if NULL, it is automatically determined.
  remove_axis_labels = FALSE, # Whether to remove the axis labels.
  colors = c("blue", "white", "red"), # Colors for the heatmap.
  # Whether to use pheatmap instead of ggplot2.
  # Some aesthetic options above do not apply when using pheatmap.
  use_pheatmap = TRUE,
  overwrite = FALSE
) {
  file_out <- fn.kmer_heatmap(
    sample_type = sample_type,
    chrom = chrom,
    strand = strand,
    ribo_nuc = ribo_nuc,
    window_side = window_side,
    window_radius = window_radius,
    kmer_size = kmer_size,
    normalize_method = normalize_method
  )

  if (!overwrite && file.exists(file_out)) {
    message("File already exists: ", file_out)
    return(invisible())
  }

  sample_list <- const.get_sample_list(sample_type)
  kmer_list <- kmer.get_all_kmers(kmer_size)

  data <- (
    kmer_heatmaps.load_data(
      sample_list = sample_list,
      chrom = chrom,
      strand = strand,
      ribo_nuc = ribo_nuc,
      window_side = window_side,
      window_radius = window_radius,
      kmer_size = kmer_size
    ) |>
    dplyr::mutate(
      kmer = factor(
        kmer,
        if (rev_kmer_list) rev(kmer_list) else kmer_list,
        ordered = TRUE
      )
    )
  )

  val_range <- if (normalize_method == "ratio") {
    c(0, 2)
  } else {
    c(0, 2 / length(kmer_list))
  }

  if (use_pheatmap) {
    library(pheatmap)
    data <- (
      data |>
      tidyr::pivot_wider(
        names_from = sample,
        values_from = value,
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
      width = 1 + (ncol(mat) * 0.3),
      height = 1 + (nrow(mat) * 0.35),
      unit = "in",
      res = 600
    )
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      breaks = seq(val_range[1], val_range[2], length.out = 200),
      fontsize_number = 10,
      number_color = "black",
      fontsize_row = 12,
      fontsize_col = 12,
      color = grDevices::colorRampPalette(colors)(200),
      labels_row = rownames(mat),
      labels_col = colnames(mat)
    )
    dev.off()
  } else {
    data <- (
      data |>
      dplyr::mutate(
        sample = const.get_sample_label_factor(sample, sample_list)
      )
    )

    if (is.null(aspect_ratio)) {
      aspect_ratio <- length(kmer_list) / length(sample_list)
    }

    ticks <- c(
      val_range[[1]],
      (val_ranges[[1]] + val_ranges[[2]]) / 2,
      val_range[[2]]
    )
    ticks_labels <- sprintf("%g", ticks)
    ticks_colors <- colors
    width <- max(length(sample_list), 8)
    height <- max(length(kmer_list), 8)

    margin <- c(
      t = 3, # inches
      r = 3, # inches
      b = 0, # inches
      l = 0 # inches
    )

    (
      data |>
      ggplot2::ggplot(ggplot2::aes(x = sample, y = kmer)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value)) +
      ggplot2::geom_text(
        ggplot2::aes(label = format(value, digits = 2, nsmall = 2)),
        size = freq_font_size,
        hjust = "middle",
        vjust = "middle"
      ) +
      ggplot2::labs(
        title = paste0(
          "Kmer heatmap\n",
          "Sample type = ", sample_type, "\n",
          "Chrom = ", chrom, "\n",
          "Strand = ", strand, "\n",
          "Ribo nucs = ", paste0(ribo_nuc, collapse = " "), "\n",
          "Window side = ", window_side, "\n",
          "Window radius = ", utils.get_pretty_base_pairs(window_radius), "\n",
          "Kmer size = ", kmer_size, "\n",
          "Normalize method = ", normalize_method, "\n",
          "X-axis shows samples\n",
          "Y-axis shows kmers"
        )
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::scale_fill_gradientn(
        limits = c(val_range[[1]], val_range[[2]]),
        breaks = ticks,
        labels = ticks_labels,
        colors = ticks_colors,
        oob = scales::squish
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_colorbar(
          title = (
            if (normalize_method == "ratio") {
              "Frequency Ratio\n(Sample / Background)"
            } else if(normalize_method == "sum1") {
              "Normalized Frequency\n(Sample / Background)"
            } else if (normalize_method == "none") {
              "Frequency"
            } else {
              stop("Unknown normalize method: ", normalize_method)
            }
          ),
          barheight = ggplot2::unit(3, "in"),
          barwidth = ggplot2::unit(0.5, "in"),
          ticks.colour = "black",
          ticks.linewidth = 1,
          frame.colour = "black",
          frame.linewidth = 1
        )
      ) +
      ggplot2::theme_classic(base_size = base_font_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "right",
        axis.title = if (remove_axis_labels) {
          ggplot2::element_blank()
        } else {
          ggplot2::element_text()
        },
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        aspect.ratio = aspect_ratio
      )
    ) |>
    utils.write_ggplot(
      file = file_out,
      width = width,
      height = height + settings.PLOT_SAMPLE_LABEL_MARGIN,
      limitsize = FALSE
    )
  }
}

kmer_heatmaps.do_main <- function(test = FALSE, overwrite = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    kmer_size_list <- c(1, 2)
    window_radius_list <- 50
    window_side_list <- "both"
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ_DRAW
    )
    kmer_size_list <- c(1, 2)
    window_radius_list <- 50
    window_side_list <- "both"
  }

  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      chrom = const.CHROM_ALL,
      strand = const.STRAND_BOTH,
      ribo_nuc = const.RIBO_NUC_ALL,
      window_side = window_side_list,
      window_radius = window_radius_list,
      kmer_size = kmer_size_list,
      normalize_method = "ratio"
    ),
    function(
      sample_type,
      chrom,
      strand,
      ribo_nuc,
      window_side,
      window_radius,
      kmer_size,
      normalize_method
    ) {
      kmer_heatmaps.plot(
        sample_type = sample_type,
        chrom = chrom,
        strand = strand,
        ribo_nuc = ribo_nuc,
        window_side = window_side,
        window_radius = window_radius,
        kmer_size = kmer_size,
        normalize_method = normalize_method,
        base_font_size = settings.BASE_FONT_SIZE,
        overwrite = overwrite
      )
    }
  )
}

kmer_heatmaps.do_precompute <- function(test = FALSE, overwrite = FALSE) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    kmer_size_list <- c(1, 2)
    window_radius_list <- 50
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ_DRAW
    )
    kmer_size_list <- c(1, 2)
    window_radius_list <- 50
  }
  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      window_radius = window_radius_list,
      kmer_size = kmer_size_list
    ),
    function(sample_type, window_radius, kmer_size) {
      kmer_heatmaps.save_kmer_data(
        sample_list = const.get_sample_list(sample_type),
        window_radius = window_radius,
        kmer_size = kmer_size,
        overwrite = overwrite
      )
    }
  )
  for (kmer_size in kmer_size_list) {
    kmer_heatmaps.save_kmer_data_bg(
      kmer_size = kmer_size,
      overwrite = overwrite
    )
  }
}
