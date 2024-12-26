# Code for making the publication figures.
# This will create a subset of the figures in "run_analysis.R"
# and write them to a directory suffixed with "_publication"

# Initialize
source("channagiri/init.R")
# Add suffix to the output path.
init.init(
  no_title = TRUE,
  plot_extra_suffix = "_publication",
  plot_dpi = 400,
  plot_sample_label_margin = 10
)

# Source the analysis functions
source("channagiri/analysis_summary.R")
source("channagiri/analysis_kmer_heatmaps.R")
source("channagiri/analysis_kmer_correlation.R")
source("channagiri/analysis_kmer_hist_offset.R")
source("channagiri/analysis_distribution.R")
source("channagiri/analysis_line_graph.R")
source("channagiri/analysis_neighborhood_counts.R")
source("channagiri/analysis_between_ribo_dist.R")
source("channagiri/analysis_dnaseq.R")
source("channagiri/analysis_overlap_clusters.R")

# Distribution analysis figures
figures.do_distribution <- function() {
  sample_type_list_list <- list(
    c(
      "Sample" = const.SAMPLE_NORMAL,
      "Random" = const.SAMPLE_UNIFORM
    ),
    c(
      "Sample" = const.SAMPLE_NORMAL,
      "Random" = const.SAMPLE_UNIFORM_CHRALL
    )
  )

  tidyr::expand_grid(
    chrom = c(const.CHROMS, const.CHROM_ALL),
    sample_type = sample_type_list_list,
    window_width = 1e5
  ) |>
  purrr::pwalk(
    function(chrom, sample_type_list, window_width) {
      distribution.plot_density(
        sample_type_list = unname(sample_type_list),
        chrom = chrom,
        strand = const.STRAND_BOTH,
        ribo_nuc = const.RIBO_NUC_ALL,
        window_width = window_width,
        normalize_method = "none",
        remove_gaps = TRUE,
        p_cutoff = 0.99,
        base_font_size = 60,
        sample_type_name_list = names(sample_type_list),
        remove_legend = TRUE,
        remove_axis_labels = TRUE,
        n_ticks_y = 3
      )
      distribution.make_test_table(
        sample_type_list = unname(sample_type_list),
        chrom = chrom,
        strand = const.STRAND_BOTH,
        ribo_nuc = const.RIBO_NUC_ALL,
        window_width = window_width,
        normalize_method = "none",
        remove_gaps = TRUE,
        exact = FALSE,
        alternative = "two.sided"
      )
    }
  )
}

figures.do_kmer_correlation <- function() {
  sample_type_list <- c(
    const.SAMPLE_NORMAL,
    const.SAMPLE_UNIFORM,
    const.SAMPLE_UNIFORM_CHRALL
  )

  tidyr::expand_grid(
    sample_type = sample_type_list,
    chrom = c(const.CHROMS, const.CHROM_ALL),
    strand = const.STRAND_BOTH,
    ribo_nuc = const.RIBO_NUC_ALL
  ) |>
  purrr::pwalk(
    function(sample_type, chrom, strand, ribo_nuc) {
      tidyr::expand_grid(
        window_width = c(1e3, 1e4, 1e5, 1e6),
        kmer_size = c(1, 2)
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
            base_font_size = 30,
            cor_font_size = 8,
            rev_kmer_list = TRUE,
            aspect_ratio = if (kmer_size == 1) {
              1 / 8
            } else if (kmer_size == 2) {
              1 / 2
            } else {
              stop("Invalid kmer_size: ", kmer_size)
            },
            remove_axis_labels = TRUE,
            use_pheatmap = TRUE
          )
        }
      )
    }
  )
}

figures.do_kmer_heatmaps <- function(overwrite = FALSE) {
  sample_type_list <- c(
    const.SAMPLE_NORMAL,
    const.SAMPLE_UNIFORM,
    const.SAMPLE_UNIFORM_CHRALL
  )

  purrr::pwalk(
    tidyr::expand_grid(
      sample_type = sample_type_list,
      chrom = c(const.CHROMS, const.CHROM_ALL),
      strand = const.STRAND_BOTH,
      ribo_nuc = const.RIBO_NUC_ALL,
      window_side = "both",
      window_radius = 50,
      kmer_size = c(1, 2),
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
        base_font_size = 30,
        freq_font_size = 8,
        rev_kmer_list = TRUE,
        aspect_ratio = if (kmer_size == 1) {
          1 / 8
        } else if (kmer_size == 2) {
          1 / 2
        } else {
          stop("Invalid kmer_size: ", kmer_size)
        },
        remove_axis_labels = TRUE,
        overwrite = overwrite
      )
    }
  )
}

figures.do_main <- function() {
  figures.do_distribution()
  figures.do_kmer_correlation()
  figures.do_kmer_heatmaps(overwrite = TRUE)
}

figures.do_main()
