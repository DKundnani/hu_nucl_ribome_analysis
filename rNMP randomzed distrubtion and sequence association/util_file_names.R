# Utility functions for generating files names

fn.NO_REPEAT_SUFFIX <- if (settings.NO_REPEAT) "_no_repeat" else ""
fn.NO_TITLE_SUFFIX <- if (settings.NO_TITLE) "_no_title" else ""
fn.OUTPUT_EXTRA_SUFFIX <- settings.OUTPUT_EXTRA_SUFFIX
fn.PLOT_EXTRA_SUFFIX <- settings.PLOT_EXTRA_SUFFIX
fn.PLOT_FILE_TYPE_SUFFIX <- paste0("_", settings.PLOT_FILE_TYPE)

fn.DATA_DIR <- file.path(
  settings.DATA_DIR_BASE,
  "channagiri",
  paste0("data", fn.NO_REPEAT_SUFFIX)
)
fn.DATA_COMMON_DIR <- file.path(
  settings.DATA_DIR_BASE,
  "channagiri",
  "data_common"
)
fn.PLOTS_DIR <- file.path(
  settings.DATA_DIR_BASE,
  "channagiri",
  paste0(
    "plot",
    fn.NO_REPEAT_SUFFIX,
    fn.NO_TITLE_SUFFIX,
    fn.PLOT_EXTRA_SUFFIX,
    fn.PLOT_FILE_TYPE_SUFFIX
  )
)
fn.DATA_OUTPUT_DIR <- file.path(
  settings.DATA_DIR_BASE,
  "channagiri",
  paste0(
    "data_output",
    fn.NO_REPEAT_SUFFIX,
    fn.OUTPUT_EXTRA_SUFFIX
  )
)

fn.HG38_GAPS <- file.path(fn.DATA_OUTPUT_DIR, "hg38_annotations", "hg38_gaps.rds")
fn.HG38_REPEATS <- file.path(fn.DATA_OUTPUT_DIR, "hg38_annotations", "hg38_repeats.rds")

fn.make_file_name <- function(dir, ext = NULL, ...) {
  if (is.list(dir)) {
    dir <- (
      dir |>
      purrr::discard(is.null) |>
      do.call(file.path, args = _)
    )
  }
  args_list <- list(...)
  args_list <- args_list |> purrr::discard(is.null)
  name <- stringr::str_c(names(args_list), args_list, collapse = "_")
  if (!is.null(ext)) {
    ext_suffix <- stringr::str_c(".", ext)
  } else {
    ext_suffix <- ""
  }
  file.path(dir, stringr::str_c(name, ext_suffix))
}

fn.ribos_tsv <- function(sample) {
  fn.make_file_name(
    dir = file.path(fn.DATA_DIR, "ribos"),
    sample,
    ext = "tsv"
  )
}

fn.ribos_rds <- function(sample, chrom) {
  fn.make_file_name(
    dir = file.path(fn.DATA_OUTPUT_DIR, "ribos", as.character(sample)),
    as.character(sample),
    as.character(chrom),
    ext = "rds"
  )
}

fn.kmer_heatmap <- function(
  sample_type,
  chrom,
  strand,
  ribo_nuc,
  window_side,
  window_radius,
  kmer_size,
  normalize_method
) {
  fn.make_file_name(
    dir = file.path(
      fn.PLOTS_DIR,
      "kmer_heatmap",
      as.integer(kmer_size),
      sample_type
    ),
    sample_type,
    chrom,
    strand,
    ribo_nuc,
    window_side,
    as.integer(window_radius),
    as.integer(kmer_size),
    normalize_method,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.kmer_heatmap_data <- function(
  sample,
  chrom,
  window_radius,
  kmer_size
) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "kmer_heatmap",
      as.integer(kmer_size),
      as.character(sample)
    ),
    as.character(sample),
    as.character(chrom),
    as.integer(window_radius),
    as.integer(kmer_size),
    ext = "rds"
  )
}

fn.kmer_heatmap_data_bg <- function(chrom, kmer_size) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "kmer_heatmap",
      as.integer(kmer_size),
      "BG"
    ),
    "BG",
    as.character(chrom),
    as.integer(kmer_size),
    ext = "rds"
  )
}

fn.kmer_hist_offset_plot <- function(
  sample_type,
  chrom,
  strand,
  window_size,
  kmer_size,
  normalize_method
) {
  fn.make_file_name(
    dir = file.path(
      fn.PLOTS_DIR,
      "kmer_hist_offset",
      as.integer(kmer_size)
    ),
    sample_type,
    chrom,
    strand,
    as.integer(window_size),
    as.integer(kmer_size),
    normalize_method,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.neighborhood_counts_plot <- function(
  sample_type,
  chrom,
  strand,
  neighborhood_size,
  count_type
) {
  fn.make_file_name(
    dir = file.path(
      fn.PLOTS_DIR,
      "neighborhood_counts",
      count_type,
      sample_type
    ),
    sample_type,
    as.character(chrom),
    as.character(strand),
    as.integer(neighborhood_size),
    count_type,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.neighborhood_counts_table <- function(
  sample_type,
  chrom,
  strand,
  neighborhood_size,
  count_type
) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "neighborhood_counts",
      count_type
    ),
    sample_type,
    as.character(chrom),
    as.character(strand),
    as.integer(neighborhood_size),
    ext = "tex"
  )
}

fn.between_ribo_dist_plot <- function(sample_type, chrom, strand, prob_type) {
  fn.make_file_name(
    dir = file.path(fn.PLOTS_DIR, "between_ribo_dist", sample_type),
    sample_type,
    chrom,
    strand,
    prob_type,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.overlap_clusters_plot <- function(sample_type, chrom, strand, col_name) {
  fn.make_file_name(
    file.path(fn.PLOTS_DIR, "overlap_clusters", sample_type),
    sample_type,
    chrom,
    strand,
    col_name,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.distribution_table <- function(
  prefix,
  sample_type_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  remove_gaps,
  ext,
  ...
) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "distribution_analysis",
      utils.get_pretty_base_pairs(window_width),
      paste0(sample_type_list, collapse = "_"),
      chrom
    ),
    sample_type_list[[1]],
    sample_type_list[[2]],
    prefix,
    chrom,
    strand,
    ribo_nuc,
    utils.get_pretty_base_pairs(window_width),
    normalize_method,
    if (remove_gaps) "noGaps" else "yesGaps",
    ...,
    ext = ext
  )
}

fn.distribution_plot <- function(
  prefix,
  sample_type_list,
  chrom,
  strand,
  ribo_nuc,
  window_width,
  normalize_method,
  remove_gaps,
  p_cutoff,
  sample,
  ...
) {
  fn.make_file_name(
    dir = file.path(
      fn.PLOTS_DIR,
      "distribution_analysis",
      utils.get_pretty_base_pairs(window_width),
      paste0(sample_type_list, collapse = "_"),
      chrom
    ),
    sample_type_list[[1]],
    sample_type_list[[2]],
    prefix,
    chrom,
    strand,
    ribo_nuc,
    utils.get_pretty_base_pairs(window_width),
    normalize_method,
    if (remove_gaps) "noGaps" else "yesGaps",
    p_cutoff,
    sample,
    ...,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.line_graph_plot <- function(sample_type, window_width, dir, prefix, ...) {
  fn.make_file_name(
    dir = file.path(
      fn.PLOTS_DIR,
      "line_graph",
      utils.get_pretty_base_pairs(window_width),
      dir
    ),
    prefix,
    sample_type,
    utils.get_pretty_base_pairs(window_width),
    ...,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.summary <- function(ext, ...) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "summary"
    ),
    ext = ext,
    ...
  )
}

fn.kmer_correlation_plot <- function(window_width, sample_type, prefix, kmer_size, ...) {
  fn.make_file_name(
    dir = file.path(
      fn.PLOTS_DIR,
      "kmer_correlation",
      utils.get_pretty_base_pairs(window_width),
      prefix,
      sample_type,
      as.integer(kmer_size)
    ),
    sample_type,
    prefix,
    utils.get_pretty_base_pairs(window_width),
    as.integer(kmer_size),
    ...,
    ext = settings.PLOT_FILE_TYPE
  )
}

fn.dnaseq_bed <- function(enzyme, assembly, chrom) {
  if (enzyme == "F") enzyme <- "dsF"
  file.path(
    fn.DATA_DIR,
    "dnaseq",
    paste0(
      "cd4t_",
      enzyme,
      "_",
      assembly,
      "_coverage"
    ),
    paste0(chrom, ".bed")
  )
}

fn.dnaseq_analysis_helper <- function(
  dir_main,
  ext,
  sample_type_ribo = NULL,
  prefix,
  assembly,
  chrom,
  strand = NULL,
  window_width = NULL,
  ...,
  normalize_method,
  remove_outliers
) {
  bp_str <- if (is.null(window_width)) {
    NULL
  } else  {
    utils.get_pretty_base_pairs(window_width)
  }
  fn.make_file_name(
    dir = list(
      dir_main,
      "dnaseq_analysis",
      bp_str,
      prefix,
      sample_type_ribo
    ),
    ext = ext,
    sample_type_ribo,
    prefix,
    assembly,
    chrom,
    strand,
    bp_str,
    ...,
    normalize_method,
    if (remove_outliers) "noOutliers" else "outliers"
  )
}

fn.dnaseq_analysis_plot <- function(
  sample_type_ribo = NULL,
  prefix,
  assembly,
  chrom,
  strand = NULL,
  window_width = NULL,
  ...,
  normalize_method,
  remove_outliers
) {
  fn.dnaseq_analysis_helper(
    dir_main = fn.PLOTS_DIR,
    ext = settings.PLOT_FILE_TYPE,
    sample_type_ribo = sample_type_ribo,
    prefix = prefix,
    assembly = assembly,
    chrom = chrom,
    strand = strand,
    window_width = window_width,
    ...,
    normalize_method = normalize_method,
    remove_outliers = remove_outliers
  )
}

fn.dnaseq_analysis_data <- function(
  ext,
  sample_type_ribo = NULL,
  prefix,
  assembly,
  chrom,
  strand = NULL,
  window_width = NULL,
  ...,
  normalize_method = NULL,
  remove_outliers
) {
  fn.dnaseq_analysis_helper(
    dir_main = fn.DATA_OUTPUT_DIR,
    ext = ext,
    sample_type_ribo = sample_type_ribo,
    prefix = prefix,
    assembly = assembly,
    chrom = chrom,
    strand = strand,
    window_width = window_width,
    ...,
    normalize_method = normalize_method,
    remove_outliers = remove_outliers
  )
}

fn.kmer_binned <- function(chrom, window_width, kmer_size) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "kmer_binned"
    ),
    ext = "rds",
    as.character(chrom),
    utils.get_pretty_base_pairs(window_width),
    as.integer(kmer_size)
  )
}

fn.weights_binned <- function(chrom, window_width) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "weights_binned",
      utils.get_pretty_base_pairs(window_width)
    ),
    ext = "rds",
    as.character(chrom),
    utils.get_pretty_base_pairs(window_width)
  )
}

fn.ribos_binned <- function(sample, chrom, window_width) {
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "ribos_binned",
      utils.get_pretty_base_pairs(window_width),
      as.character(sample)
    ),
    ext = "rds",
    as.character(sample),
    as.character(chrom),
    utils.get_pretty_base_pairs(window_width)
  )
}

fn.dnaseq_rds <- function(sample, assembly, chrom) {
  if (sample == "F") sample <- "dsF"
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "dnaseq",
      paste0(
        "cd4t_",
        sample,
        "_",
        assembly
      )
    ),
    chrom,
    ext = "rds"
  )
}

fn.dnaseq_binned <- function(sample, assembly, chrom, window_width) {
  if (sample == "F") sample <- "dsF"
  fn.make_file_name(
    dir = file.path(
      fn.DATA_OUTPUT_DIR,
      "dnaseq_binned",
      utils.get_pretty_base_pairs(window_width),
      paste0(
        "cd4t_",
        sample,
        "_",
        assembly
      )
    ),
    paste0(
      "cd4t_",
      sample,
      "_",
      assembly
    ),
    chrom,
    utils.get_pretty_base_pairs(window_width),
    ext = "rds"
  )
}
