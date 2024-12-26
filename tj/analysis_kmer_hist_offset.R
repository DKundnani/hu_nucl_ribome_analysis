kmer_hist_offset.get_data <- function(
  sample_list,
  chrom,
  window_size,
  kmer_size,
  normalize_method = "ratio" # choices = "none" "ratio", "sum1"
) {
  purrr::pmap(
    tidyr::expand_grid(
      sample = sample_list,
      chrom = chrom
    ),
    function(sample, chrom) {
      utils.load_ribos_rds(sample = sample, chrom = chrom) |>
      x => list( # Split by strand
        strand = c("+", "-"),
        granges = list(
          x[GenomicRanges::strand(x) == "+",],
          x[GenomicRanges::strand(x) == "-",]
        )
      ) |>
      purrr::pmap(
        function(strand, granges) {
          GenomicRanges::GRanges(
            seqnames = chrom,
            ranges = IRanges::IRanges(
              start = GenomicRanges::start(granges) - window_size,
              end = GenomicRanges::start(granges) + window_size
            ),
            strand = GenomicRanges::strand(granges)
          ) |>
          x => Biostrings::getSeq(
            BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
            x
          ) |>
          as.character() |>
          x => utils.get_kmer_hist_by_offset(
            seq = x,
            count = granges$count,
            kmer_size = kmer_size,
            zero_offset = window_size + 1
          ) |>
          dplyr::mutate(
            sample = !!sample,
            chrom = !!chrom,
            strand = !!strand,
            .before = 1
          )
        }
      ) |>
      purrr::list_rbind() |>
      tidyr::pivot_longer(
        dplyr::all_of(kmer.get_all_kmers(kmer_size)),
        names_to = "kmer",
        values_to = "freq"
      ) |>
      dplyr::group_by(
        sample,
        chrom,
        strand,
        offset
      ) |>
      dplyr::mutate(freq = freq / sum(freq)) |>
      dplyr::ungroup()
    }
  ) |>
  purrr::list_rbind()
}

kmer_hist_offset.get_data_bg <- function(chrom, window_size, kmer_size) {
  tidyr::expand_grid(
    chrom = chrom,
    strand = const.STRANDS,
    tibble::tibble(
      start = 1,
      end = utils.get_chrom_sizes(!!chrom)
    )
  ) |>
  data => dplyr::bind_cols(
    data[, c("chrom", "strand")],
    utils.get_kmer_freqs(
      data[["chrom"]],
      data[["strand"]],
      data[["start"]],
      data[["end"]],
      kmer_size
    )
  ) |>
  tidyr::expand_grid(offset = seq(-window_size, window_size + kmer_size - 1)) |>
  tidyr::pivot_longer(
    dplyr::all_of(kmer.get_all_kmers(kmer_size)),
    names_to = "kmer",
    values_to = "freq_bg"
  )
}

kmer_hist_offset.plot <- function(
  sample_type,
  chrom,
  window_size,
  kmer_size,
  normalize_method = "ratio", # choices = "none", "ratio" ,"sum1"
  base_font_size = 11
) {
  kmer_list <- kmer.get_all_kmers(kmer_size)
  sample_list <- const.get_sample_list(sample_type)

  data <- kmer_hist_offset.get_data(
    sample_list = sample_list,
    chrom = chrom,
    window_size = window_size,
    kmer_size = kmer_size,
    normalize_method = normalize_method
  )
  if (normalize_method %in% c("ratio", "sum1")) {
    data_bg <- kmer_hist_offset.get_data_bg(
      chrom = chrom,
      window_size = window_size,
      kmer_size = kmer_size
    )
    data <- (
      dplyr::left_join(
        data,
        data_bg,
        by = c("chrom", "strand", "offset", "kmer")
      ) |>
      dplyr::mutate(freq = freq / freq_bg) |>
      dplyr::select(!freq_bg) |>
      dplyr::group_by(sample, chrom, strand, offset) |>
      dplyr::mutate(
        freq = if (normalize_method == "sum1") {
          freq / sum(freq)
        } else {
          freq
        }
      ) |>
      dplyr::ungroup()
    )
  } else if (normalize_method == "none") {
  } else {
    stop(paste0("Unknown normalize method: ", normalize_method))
  }

  kmer_list <- kmer.get_all_kmers(kmer_size)
  offset_list <- seq(-window_size, window_size + kmer_size - 1)
  y_lim <- if (normalize_method == "ratio") {
    c(0, 2)
  } else {
    c(0, 2 / length(kmer_list))
  }
  tick_breaks <- if (normalize_method == "ratio") {
    c(0, 1, 2)
  } else {
    c(0, 1 / length(kmer_list), 2 / length(kmer_list))
  }
  tick_labels <- if (normalize_method == "ratio") {
    c("0", "1", "2")
  } else {
    format(tick_breaks, digits = 2, nsmall = 2)
  }
  hline_y <- if (normalize_method == "ratio") {
    1
  } else {
    1 / length(kmer_list)
  }

  data |>
  dplyr::nest_by(strand, .key = "data_1") |>
  purrr::pwalk(
    function(strand, data_1) {
      (
        data_1 |>
        dplyr::mutate(
          sample = const.get_sample_label_factor(sample, sample_list)
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = offset, y = freq, color = kmer)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = hline_y, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        ggplot2::coord_cartesian(ylim = y_lim) +
        ggplot2::scale_y_continuous(
          breaks = tick_breaks,
          labels = tick_labels
        ) +
        ggplot2::scale_x_continuous(
          breaks = c(
            -window_size,
            -window_size / 2,
            0,
            window_size / 2,
            window_size
          ),
          labels = c(
            -window_size,
            -window_size / 2,
            0,
            window_size / 2,
            window_size
          )
        ) +
        ggplot2::labs(
          title = paste0(
            "Kmer offset histogram\n",
            "Sample type = ", sample_type, "\n",
            "Chromosome = ", chrom, "\n",
            "Strand = ", strand, "\n",
            "Window size = ", utils.get_pretty_base_pairs(window_size), "\n",
            "Kmer size = ", kmer_size, "\n",
            "Normalize method = ", normalize_method, "\n",
            "X-axis shows offset from ribo\n",
            "Y-axis shows (normalized) frequency of kmer at offset"
          )
        ) +
        ggplot2::xlab("Offset from ribo") +
        ggplot2::ylab(
          if (normalize_method == "sum1") {
            "Normalized frequency"
          } else if (normalize_method == "ratio") {
            "Frequency ratio (freq / freq_bg)"
          } else {
            "Frequency"
          }
        ) +
        ggplot2::theme_classic(base_size = base_font_size) +
        ggplot2::theme(
          plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          strip.text = element_text(size = 7),
          strip.background = element_rect(color = "black"),
          legend.key.size = unit(1, "inches")
        ) +
        ggplot2::facet_wrap(dplyr::vars(sample))
      ) |>
      utils.write_ggplot(
        fn.kmer_hist_offset_plot(
          sample_type = sample_type,
          chrom = chrom,
          strand = strand,
          window_size = window_size,
          kmer_size = kmer_size,
          normalize_method = normalize_method
        ),
        width = max(ceiling(sqrt(length(sample_list))) * 3, 4),
        height = max(ceiling(sqrt(length(sample_list))) * 2, 4),
        limitsize = FALSE
      )
    }
  )
}

kmer_hist_offset.do_main <- function(
  test = FALSE,
  normalize_method = "ratio"
) {
  if (test) {
    sample_type_list <- const.SAMPLE_TEST
    chrom_list <- "chr1"
  } else {
    sample_type_list <- c(
      const.SAMPLE_NORMAL,
      const.SAMPLE_SHUFFLE,
      const.SAMPLE_UNIFORM,
      const.SAMPLE_UNIFORM_CHRALL,
      const.SAMPLE_DNASEQ_DRAW
    )
    chrom_list <- const.CHROMS
  }

  tidyr::expand_grid(
    sample_type = sample_type_list,
    chrom = chrom_list,
    window_size = 50L,
    kmer_size = 1L,
    normalize_method = normalize_method
  ) |>
  purrr::pwalk(
    function(sample_type, chrom, window_size, kmer_size, normalize_method) {
      kmer_hist_offset.plot(
        sample_type = sample_type,
        chrom = chrom,
        window_size = window_size,
        kmer_size = kmer_size,
        normalize_method = normalize_method,
        base_font_size = settings.BASE_FONT_SIZE
      )
    }
  )
}
