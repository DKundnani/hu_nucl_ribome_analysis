# common utility functions

utils.create_parent_dir <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

utils.log <- function(...) {
  if (settings.LOG_ENABLE) {
    args <- list(...)
    args_names <- names(args)
    if (is.null(args_names)) {
      args_names <- rep(NA, length(args))
    }
    arg_names <- ifelse(is.na(args_names), "", stringr::str_c(args_names, "="))
    args <- stringr::str_c(stringr::str_c(arg_names, args), collapse = " ")
    cat(
      stringr::str_c(
        stringr::str_split(lubridate::now(), " ")[[1]][2],
        ": ",
        args,
        "\n"
      )
    )
  }
}

utils.enable_log <- function() {
  settings.LOG_ENABLE <<- TRUE
}

utils.disable_log <- function() {
  settings.LOG_ENABLE <<- FALSE
}

# Ignore logging for the evaluation of the expression.
# Note: forces evaulation of the expression.
utils.ignore_log <- function(x) {
  settings.LOG_ENABLE <<- FALSE
  force(x)
  settings.LOG_ENABLE <<- TRUE
  x
}

utils.read_delim <- function(file, delim, ...) {
  utils.log("(input)", file)
  readr::read_delim(
    file,
    delim = delim,
    show_col_types = FALSE,
    progress = FALSE,
    ...
  )
}

utils.read_tsv <- function(file, ...) {
  utils.read_delim(file, "\t", ...)
}

utils.read_csv <- function(file, ...) {
  utils.read_delim(file, ",", ...)
}

utils.read_rds <- function(file, ...) {
  utils.log("(input)", file)
  readRDS(file, ...)
}

utils.write_delim <- function(data, file, delim, na = "", ...) {
  force(data)
  utils.log("(output)", file)
  utils.create_parent_dir(file)
  readr::write_delim(
    data,
    file,
    delim = delim,
    na = na,
    progress = FALSE,
    ...
  )
}

utils.write_tsv <- function(data, file, ...) {
  utils.create_parent_dir(file)
  utils.write_delim(data, file, "\t", ...)
}

utils.write_csv <- function(data, file, ...) {
  utils.create_parent_dir(file)
  utils.write_delim(data, file, ",", ...)
}

utils.write_bedgraph <- function(grange, file, ...) {
  force(grange)
  utils.log("(output)", file)
  utils.create_parent_dir(file)
  rtracklayer::export.bedGraph(
    grange,
    file,
    ...
  )
}

utils.write_bed <- function(grange, file, ...) {
  force(grange)
  utils.log("(output)", file)
  utils.create_parent_dir(file)
  rtracklayer::export.bed(
    grange,
    file,
    ...
  )
}

utils.write_rds <- function(data, file, ...) {
  force(data)
  utils.log("(output)", file)
  utils.create_parent_dir(file)
  saveRDS(
    data,
    file,
    compress = TRUE,
    ...
  )
}

# BED formats have 0-based coordinates, with the start coordinate
# inclusive and the end coordinate exclusive.
# Convert this to 1-based coordinates, with both coordinates inclusive.
utils.format_bed_as_table <- function(bed_data) {
  # no change needed for end coordinates
  bed_data |> dplyr::mutate(start = start + 1)
}

utils.read_bed <- function(file, ...) {
  utils.log("(input)", file)

  rtracklayer::import.bed(file, ...) |>
  tibble::as_tibble() |>
  dplyr::rename(chrom = seqnames) |>
  dplyr::mutate(
    chrom = as.character(chrom),
    strand = as.character(strand)
  )
}

utils.read_bedgraph <- function(file, ...) {
  utils.log("(input)", file)

  rtracklayer::import.bedGraph(file, ...) |>
  tibble::as_tibble() |>
  dplyr::rename(chrom = seqnames) |>
  dplyr::mutate(strand = as.character(strand))
}

utils.write_ggplot <- function(
  plot,
  file,
  height,
  width,
  margin = NULL,
  ...
) {
  force(plot)
  utils.log("(output)", file)
  utils.create_parent_dir(file)
  if (settings.NO_TITLE) {
    plot <- plot + ggplot2::theme(plot.title = ggplot2::element_blank())
  }
  ggplot2::ggsave(
    file,
    plot = plot,
    dpi = settings.PLOT_DPI,
    height = height + if (!settings.NO_TITLE) settings.PLOT_TITLE_MARGIN else 0,
    width = width,
    units = "in",
    ...
  )
}

utils.write_xtable <- function(table, file, ...) {
  force(table)
  utils.log("(output)", file)
  utils.create_parent_dir(file)
  print(
    table,
    file = file,
    include.rownames = FALSE,
    comment = FALSE,
    timestamp = FALSE,
    ...
  )
}

utils.get_chrom_sizes <- function(chroms = const.CHROMS) {
  GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[chroms]
}

utils.get_hg38_gaps <- function() {
  utils.read_rds(fn.HG38_GAPS)
}

utils.granges_to_tibble <- function(x) {
  x |>
  tibble::as_tibble() |>
  dplyr::rename(chrom = seqnames) |>
  dplyr::mutate(
    chrom = as.character(chrom),
    strand = as.character(strand)
  )
}

utils.load_ribos_rds <- function(sample, chrom) {
  chrom_list <- const.expand_chrom(chrom)

  purrr::map(
    chrom_list,
    function(chrom) {
      fn.ribos_rds(sample = sample, chrom = chrom) |>
      utils.read_rds()
    }
  ) |>
  x => do.call(plyranges::bind_ranges, x)
}

utils.load_dnaseq_rds <- function(sample, chrom, assembly = "hg38") {
  chrom_list <- const.expand_chrom(chrom)

  purrr::map(
    chrom_list,
    function(chrom) {
      fn.dnaseq_rds(
        sample = sample,
        chrom = chrom,
        assembly = assembly
      ) |>
      utils.read_rds()
    }
  ) |>
  x => do.call(plyranges::bind_ranges, x)
}

utils.load_ribos_binned <- function(
  sample,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method = "divideNonGap",
  remove_gaps = FALSE # Remove windows that are totally contained in a gap
) {
  chrom_list <- const.expand_chrom(chrom)
  strand_list <- const.expand_strand(strand)
  ribo_nuc_list <- const.expand_ribo_nuc(ribo_nuc)
  
  ranges <- (
    purrr::map(
      chrom_list,
      function(chrom) {
        ranges <- (
          fn.ribos_binned(
            sample = sample,
            chrom = chrom,
            window_width = window_width
          ) |>
          utils.read_rds()
        )
        if (is.null(ranges)) {
          ranges <- utils.get_ranges_binned(
            chrom = chrom,
            window_width = window_width,
            strand = TRUE
          )
          ranges$count = matrix(
            0,
            nrow = NROW(ranges),
            ncol = length(const.RIBO_NUCS),
            dimnames = list(NULL, const.RIBO_NUCS)
          )
        }
        ranges <- ranges[GenomicRanges::strand(ranges) %in% strand_list, ]
        ranges$count <- ranges$count[, ribo_nuc_list, drop = FALSE]
        ranges$count <- rowSums(ranges$count)

        # These two variables are only needed if
        # using "divideNonGap" normalize method or removing gaps.
        ranges_weights <- NULL
        window_index <- NULL
        if ((normalize_method == "divideNonGap") || remove_gaps) {
          ranges_weights <- utils.load_weights_binned(
            chrom = chrom,
            window_width = window_width
          )
          window_index <- (
            GenomicRanges::findOverlaps(
              ranges,
              ranges_weights,
              type = "equal",
              select = "first"
            ) |>
            as.integer()
          )
          if (any(is.na(window_index))) {
            stop("Some ranges not found in weights:", sample, chrom, window_width)
          }
        }
      
        if (normalize_method == "divideNonGap") {
          ranges$count <- ranges$count / ranges_weights$count[window_index]
          ranges$count[!is.finite(ranges$count)] <- NA
        } else if (normalize_method == "width") {
          ranges$count <- ranges$count / window_width
        } else if (normalize_method == "none") {
        } else {
          stop("Unknown normalize_method:", normalize_method)
        }

        if (remove_gaps) {
          ranges <- ranges[ranges_weights$count[window_index] > 0]
        }

        ranges
      }
    )
  ) |>
  x => do.call(plyranges::bind_ranges, x)
}

utils.load_dnaseq_binned <- function(
  sample,
  chrom,
  strand,
  window_width,
  assembly = "hg38",
  normalize_method = "divideNonGap",
  remove_outliers = FALSE
) {
  chrom_list <- const.expand_chrom(chrom)
  strand_list <- const.expand_strand(strand)
  
  ranges <- (
    purrr::map(
      chrom_list,
      function(chrom, strand) {
        ranges <- (
          fn.dnaseq_binned(
            sample = sample,
            chrom = chrom,
            window_width = window_width,
            assembly = "hg38"
          ) |>
          utils.read_rds()
        )
        ranges <- ranges[GenomicRanges::strand(ranges) %in% strand_list, ]
        if (normalize_method == "divideNonGap") {
          ranges_weights <- utils.load_weights_binned(
            chrom = chrom,
            window_width = window_width
          )
          window_index <- (
            GenomicRanges::findOverlaps(
              ranges,
              ranges_weights,
              type = "equal",
              select = "first"
            ) |>
            as.integer()
          )
          if (any(is.na(window_index))) {
            stop("Some ranges not found in weights:", sample, chrom, window_width)
          }
          ranges$count <- ranges$count / ranges_weights$count[window_index]
          ranges$count[!is.finite(ranges$count)] <- NA
        } else if (normalize_method == "width") {
          ranges$count <- ranges$count / window_width
        } else {
          stop("Unknown normalize_method:", normalize_method)
        }
        ranges
      }
    ) |>
    x => do.call(plyranges::bind_ranges, x)
  )

  ranges$outlier <- (
    ranges$count >
    (
      median(ranges$count, na.rm = TRUE) +
      (4 * mad(ranges$count, na.rm = TRUE))
    )
  )
  if (remove_outliers) {
    ranges$count[ranges$outlier] <- NA
  }

  ranges
}

utils.load_weights_binned <- function(chrom, window_width) {
  fn.weights_binned(chrom = chrom, window_width = window_width) |>
  utils.read_rds()
}

utils.get_kmer_hist_by_offset <- function(seq, count, kmer_size, zero_offset) {
  all_kmers <- kmer.get_all_kmers(kmer_size)
  seq_size <- nchar(seq[[1]])
  start_index <- seq_len(seq_size - kmer_size + 1)
  end_index <- start_index + kmer_size - 1
  seq <- sapply(
    seq,
    function(str) {
      substring(
        str,
        first = start_index,
        last = end_index
      )
    },
    USE.NAMES = FALSE
  )
  seq <- t(seq)
  offset_names <- as.character(start_index - zero_offset)
  colnames(seq) <- offset_names
  seq <- tibble::as_tibble(seq)
  seq[["count"]] <- count
  seq <- (
    tidyr::pivot_longer(
      seq,
      all_of(offset_names),
      names_to = "offset",
      values_to = "kmer"
    ) |>
    dplyr::group_by(kmer, offset) |>
    dplyr::summarize(count = sum(count), .groups = "drop")
  )
  seq <- tidyr::pivot_wider(
    seq,
    id_cols = "offset",
    names_from = "kmer",
    values_from = "count",
    values_fill = 0
  )
  missing_kmers <- all_kmers[!(all_kmers %in% colnames(seq))]
  seq[, missing_kmers] <- 0
  seq <- (
    seq |>
    dplyr::select(offset, all_of(all_kmers)) |>
    dplyr::mutate(offset = as.integer(offset)) |>
    dplyr::arrange(offset)
  )
  seq
}

utils.coverage_intervals <- function(
  chrom,
  start,
  end,
  count = NULL,
  only_orig = FALSE
) {
  if (length(start) == 0) {
    tibble::tibble(
      chrom = character(0),
      start = integer(0),
      end = integer(0),
      count = integer(0)
    )
  } else {
    if (is.null(count)) {
      count <- rep(1, length(start))
    }
    ribos_granges <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges = IRanges::IRanges(
        start = start,
        end = end,
      ),
      seqlengths = utils.get_chrom_sizes()
    )
    coverage_out <- GenomicRanges::coverage(ribos_granges, weight = count)
    counts_rle <- purrr::map_dfr(
      chrom,
      function(chrom) {
        counts_rle <- coverage_out[[chrom]]
        coverage_end <- cumsum(S4Vectors::runLength(counts_rle))
        coverage_start <- c(1, coverage_end[1 : (length(coverage_end) - 1)] + 1)
        coverage_count <- S4Vectors::runValue(counts_rle)
        if (only_orig) {
          coverage_granges <- GenomicRanges::GRanges(
            seqnames = chrom,
            ranges = IRanges::IRanges(
              start = coverage_start,
              end = coverage_end
            ),
            seqlengths = utils.get_chrom_sizes()
          )
          overlaps <- GenomicRanges::findOverlaps(
            ribos_granges,
            coverage_granges,
            select = "first"
          )
          tibble::tibble(
            chrom = chrom,
            start = start,
            end = end,
            count = coverage_count[overlaps]
          )
        } else {
          tibble::tibble(
            chrom = chrom,
            start = coverage_start,
            end = coverage_end,
            count = coverage_count
          )
        }
      }
    )
  }
}

utils.reduce_intervals <- function(start, end) {
  if (length(start) == 0) {
    tibble::tibble(start = integer(0), end = integer(0))
  } else {
    IRanges::IRanges(start = start, end = end) |>
    IRanges::reduce() |>
    x => tibble::tibble(start = start(x), end = end(x))
  }
}

utils.subtract_intervals <- function(
  start_1,
  end_1,
  start_2,
  end_2,
  total_size
) {
  ranges_1 <- IRanges::IRanges(start = start_1, end = end_1)
  ranges_2 <- IRanges::gaps(
    IRanges::IRanges(
      start = c(0, start_2, total_size + 1),
      end = c(0, end_2, total_size + 1)
    )
  )
  ranges_subtract <- IRanges::overlapsRanges(ranges_1, ranges_2)
  tibble::tibble(
    start = start(ranges_subtract),
    end = end(ranges_subtract)
  )
}

utils.intersect_intervals <- function(
  start_1,
  end_1,
  start_2,
  end_2
) {
  ranges_1 <- IRanges::IRanges(start = start_1, end = end_1)
  ranges_2 <- IRanges::IRanges(start = start_2, end = end_2)
  ranges_intersect <- IRanges::overlapsRanges(ranges_1, ranges_2)
  tibble::tibble(
    start = start(ranges_intersect),
    end = end(ranges_intersect)
  )
}

# Get empirical quantiles for given discrete distribution
utils.get_quantiles <- function(x, freq, p_list) {
  cdf <- cumsum(freq)
  cdf <- cdf / cdf[[length(cdf)]]
  purrr::map_dbl(
    p_list,
    function(p) {
      indices <- which(cdf >= p)
      if (length(indices) == 0) {
        NA
      } else {
        x[indices[[1]]]
      }
    }
  )
}

utils.get_probs <- function(x, freq, q_list) {
  purrr::map_dbl(
    q_list,
    function(q) {
      sum(freq[x <= q]) / sum(freq)
    }
  )
}

# Convenience function for debugging a pipeline
# Enters browser before passing through a value
utils.inspect <- function(x) {
  force(x)
  browser()
  x
}

utils.get_pretty_base_pairs <- function(x) {
  x |>
  purrr::map_chr(
    function(y) {
      if (y < 1e3) {
        stringr::str_c(y, "b")
      } else if (y < 1e6) {
        stringr::str_c(format(y / 1e3, digits = 1), "kb")
      } else if (y < 1e9) {
        stringr::str_c(format(y / 1e6, digits = 1), "mb")
      }
    }
  )
}

utils.get_tick_info <- function(chrom_size) {
  tick_size_list <- cumprod(rep(c(5, 2), 10))
  tick_size <- purrr::detect(tick_size_list, function(x) (chrom_size / x) < 25)
  tick_breaks <- seq(0, chrom_size %/% tick_size) * tick_size
  tick_labels <- purrr::map_chr(
    tick_breaks,
    utils.get_pretty_base_pairs
  )
  list(
    breaks = tick_breaks,
    labels = tick_labels
  )
}

utils.block_output <- function(block) {
  if (block) {
    if (.Platform$OS.type == "windows") {
      sink("NUL") # special file on Windows for no output
    } else if (.Platform$OS.type == "unix") {
      sink("/dev/null") # special file on Unix for no output
    }
  } else {
    if (.Platform$OS.type %in% c("windows", "unix")) {
      sink(NULL) # undo the sink
    }
  }
}

utils.get_kmer_counts <- function(chrom, strand, start, end, kmer_size) {
  utils.log("Getting kmer counts")
  Biostrings::oligonucleotideFrequency(
    IRanges::Views(
      BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
      GenomicRanges::GRanges(
        seqnames = chrom,
        strand = strand,
        ranges = IRanges::IRanges(start = start, end = end)
      )
    ),
    kmer_size
  )
}

utils.get_kmer_freqs <- function(chrom, strand, start, end, kmer_size) {
  x <- utils.get_kmer_counts(chrom, strand, start, end, kmer_size)
  sweep(x, 1, rowSums(x), "/")
}

# intersect ranges_1 and ranges_2 but get the metadata from ranges_1
utils.intersect_ranges <- function(ranges_1, ranges_2) {
  overlaps <- GenomicRanges::findOverlaps(ranges_1, ranges_2)
  ranges_1 <- ranges_1[queryHits(overlaps), ]
  ranges_2 <- ranges_2[subjectHits(overlaps), ]
  ranges_result <- GenomicRanges::pintersect(ranges_1, ranges_2)
  ranges_result$hits <- NULL # remove hits
  ranges_result
}

# "method" can be "rank", "mean", or "sum"
utils.normalize <- function(x, method, ...) {
  if (method == "rank") {
    rank(x, na.last = TRUE, ...)
  } else if (method == "mean") {
    x / mean(x, na.rm = TRUE, ...)
  } else if (method == "sum") {
    x / sum(x, na.rm = TRUE, ...)
  } else {
    stop(paste0("Unknown method: ", method))
  }
}

utils.get_neighborhood_counts <- function(
  sample,
  chrom,
  neighborhood_size,
  count_type,
  only_orig = TRUE
) {
  chrom_size <- utils.get_chrom_sizes(chrom)

  utils.load_ribos_rds(sample = sample, chrom = chrom) |>
  utils.granges_to_tibble() |>
  dplyr::nest_by(strand, .keep = TRUE) |>
  purrr::pmap(
    function(strand, data) {
      count <- (
        if (count_type == "multiple") {
          data[["count"]]
        } else if (count_type == "single") {
          rep(1, nrow(data))
        } else {
          stop(paste0("Unknown count type: ", count_type))
        }
      )

      data_cov <- utils.coverage_intervals(
        start = pmax(1, data[["start"]] - neighborhood_size),
        end = pmin(chrom_size, data[["start"]] + neighborhood_size),
        chrom = chrom,
        count = count,
        only_orig = only_orig
      )

      data |>
      dplyr::mutate(
        neighborhood_size = !!neighborhood_size,
        count_nb = data_cov[["count"]],
        count_type = !!count_type
      ) |>
      dplyr::select(
        sample,
        chrom,
        strand,
        start,
        end,
        ribo_nuc,
        count_type,
        neighborhood_size,
        count,
        count_nb
      )
    }
  ) |>
  purrr::list_rbind()
}

utils.get_ranges_binned <- function(chrom, window_width, strands = FALSE) {
  GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(
      start = 1,
      end = utils.get_chrom_sizes(chrom)
    )
  ) |>
  plyranges::slide_ranges(
    width = window_width,
    step = window_width
  ) |>
  GenomicRanges::granges() |>
  x => (
    if (strands) {
      c(
        x |> plyranges::mutate(strand = "+"),
        x |> plyranges::mutate(strand = "-")
      )
    } else {
      x
    }
  )
}
