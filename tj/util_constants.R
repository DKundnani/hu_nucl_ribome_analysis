# Utilities and definitions for commonly used constants.

const.CHROMS <- c(
  "chr1",
  "chr2",
  "chr3",
  "chr4",
  "chr5",
  "chr6",
  "chr7",
  "chr8",
  "chr9",
  "chr10",
  "chr11",
  "chr12",
  "chr13",
  "chr14",
  "chr15",
  "chr16",
  "chr17",
  "chr18",
  "chr19",
  "chr20",
  "chr21",
  "chr22"
  # "chrX", # omitted for now, TJ 2023-08-08
  # "chrM", # omitted for now, TJ 2023-08-08
  # "chrY", # omitted for now, TJ 2023-08-08
)
const.CHROM_ALL <- "chrAll" # all chromosomes

const.STRANDS <- c("+", "-")
const.STRAND_BOTH <- "+-" # both strands

const.RIBO_NUCS <- c("A", "C", "G", "T")
const.RIBO_NUC_ALL <- "nucAll"

const.expand_chrom <- function(chrom) {
  if (!is.null(chrom) && (chrom == const.CHROM_ALL)) {
    const.CHROMS
  } else {
    chrom
  }
}

const.expand_strand <- function(strand) {
  if (!is.null(strand) && (strand == const.STRAND_BOTH)) {
    const.STRANDS
  } else {
    strand
  }
}

const.expand_ribo_nuc <- function(ribo_nuc) {
  if (!is.null(ribo_nuc) && (ribo_nuc == const.RIBO_NUC_ALL)) {
    const.RIBO_NUCS
  } else {
    ribo_nuc
  }
}

const.SAMPLE_INFO <- list(
  FS185 = list(cell = "CD4T", enzyme = "F", cell_class = "WT"),
  FS187 = list(cell = "CD4T", enzyme = "F_RE1", cell_class = "WT"),
  FS188 = list(cell = "CD4T", enzyme = "F_RE1", cell_class = "WT"),
  FS189 = list(cell = "CD4T", enzyme = "F_RE1", cell_class = "WT"),
  FS186 = list(cell = "CD4T", enzyme = "RE1", cell_class = "WT"),
  FS197 = list(cell = "hESC-H9", enzyme = "F", cell_class = "WT"),
  FS198 = list(cell = "hESC-H9", enzyme = "F", cell_class = "WT"),
  FS201 = list(cell = "hESC-H9", enzyme = "F_RE1", cell_class = "WT"),
  FS199 = list(cell = "hESC-H9", enzyme = "RE1", cell_class = "WT"),
  FS331 = list(cell = "hESC-H9", enzyme = "RE2", cell_class = "WT"),
  FS303 = list(cell = "hESC-H9", enzyme = "RE3", cell_class = "WT"),
  FS326 = list(cell = "HEK293T-WT", enzyme = "F", cell_class = "WT"),
  FS391 = list(cell = "HEK293T-WT", enzyme = "F", cell_class = "WT"),
  FS203 = list(cell = "HEK293T-WT", enzyme = "F_RE1", cell_class = "WT"),
  FS333 = list(cell = "HEK293T-WT", enzyme = "RE2", cell_class = "WT"),
  FS305 = list(cell = "HEK293T-WT", enzyme = "RE3", cell_class = "WT"),
  FS327 = list(cell = "HEK293T-RNASEH2A-KO-T3-8", enzyme = "F", cell_class = "KO"),
  FS392 = list(cell = "HEK293T-RNASEH2A-KO-T3-8", enzyme = "F", cell_class = "KO"), 
  FS300 = list(cell = "HEK293T-RNASEH2A-KO-T3-8", enzyme = "RE1", cell_class = "KO"),
  FS309 = list(cell = "HEK293T-RNASEH2A-KO-T3-8", enzyme = "RE2", cell_class = "KO"),
  FS306 = list(cell = "HEK293T-RNASEH2A-KO-T3-8", enzyme = "RE3", cell_class = "KO"),
  FS329 = list(cell = "HEK293T-RNASEH2A-KO-T3-17", enzyme = "F", cell_class = "KO"),
  FS393 = list(cell = "HEK293T-RNASEH2A-KO-T3-17", enzyme = "F", cell_class = "KO"),
  FS301 = list(cell = "HEK293T-RNASEH2A-KO-T3-17", enzyme = "RE1", cell_class = "KO"),
  FS307 = list(cell = "HEK293T-RNASEH2A-KO-T3-17", enzyme = "RE2", cell_class = "KO"),
  FS310 = list(cell = "HEK293T-RNASEH2A-KO-T3-17", enzyme = "RE3", cell_class = "KO")
)

const.SAMPLE_INFO_DNASEQ <- list(
  F = list(cell = "CD4T", enzyme = "F", cell_class = "WT"),
  RE1 = list(cell = "CD4T", enzyme = "RE1", cell_class = "WT"),
  RE2 = list(cell = "CD4T", enzyme = "RE2", cell_class = "WT"),
  RE3 = list(cell = "CD4T", enzyme = "RE3", cell_class = "WT")
)

const.SAMPLES <- names(const.SAMPLE_INFO)
const.SAMPLES_ENZYME_F <- purrr::keep(
  const.SAMPLES,
  function(x) const.SAMPLE_INFO[[x]][["enzyme"]] == "F"
)
const.SAMPLES_WT <- purrr::keep(
  const.SAMPLES,
  function(x) const.SAMPLE_INFO[[x]][["cell_class"]] == "WT"
)
const.SAMPLES_KO <- purrr::keep(
  const.SAMPLES,
  function(x) const.SAMPLE_INFO[[x]][["cell_class"]] == "KO"
)
const.get_sample_hotspot <- function(sample) paste0(sample, "_hotspot")
const.SAMPLES_HOTSPOT <- purrr::map_chr(const.SAMPLES, const.get_sample_hotspot)

const.get_sample_shuffle <- function(sample) paste0(sample, "_shuffle")
const.SAMPLES_SHUFFLE <- purrr::map_chr(const.SAMPLES, const.get_sample_shuffle)
const.SAMPLES_SHUFFLE_ENZYME_F <- purrr::map_chr(const.SAMPLES_ENZYME_F, const.get_sample_shuffle)

const.get_sample_uniform <- function(sample) paste0(sample, "_uniform")
const.SAMPLES_UNIFORM <- purrr::map_chr(const.SAMPLES, const.get_sample_uniform)
const.SAMPLES_UNIFORM_ENZYME_F <- purrr::map_chr(const.SAMPLES_ENZYME_F, const.get_sample_uniform)

const.get_sample_uniform_chrAll <- function(sample) paste0(sample, "_uniform_chrAll")
const.SAMPLES_UNIFORM_CHRALL <- purrr::map_chr(const.SAMPLES, const.get_sample_uniform_chrAll)
const.SAMPLES_UNIFORM_CHRALL_ENZYME_F <- purrr::map_chr(const.SAMPLES_ENZYME_F, const.get_sample_uniform_chrAll)

const.get_sample_dnaseq_draw <- function(sample) paste0(sample, "_dnaseq_draw")
const.SAMPLES_DNASEQ_DRAW <- purrr::map_chr(const.SAMPLES, const.get_sample_dnaseq_draw)
const.SAMPLES_DNASEQ_DRAW_ENZYME_F <- purrr::map_chr(const.SAMPLES_ENZYME_F, const.get_sample_dnaseq_draw)

const.get_sample_test <- function(sample) paste0(sample, "_test")
const.SAMPLES_TEST <- purrr::map_chr(const.SAMPLES, const.get_sample_test)
const.SAMPLES_TEST_ENZYME_F <- purrr::map_chr(const.SAMPLES_ENZYME_F, const.get_sample_test)
const.SAMPLES_TEST_SHUFFLE <- purrr::map_chr(const.SAMPLES_TEST, const.get_sample_shuffle)
const.SAMPLES_TEST_UNIFORM <- purrr::map_chr(const.SAMPLES_TEST, const.get_sample_uniform)
const.SAMPLES_TEST_UNIFORM_CHRALL <- purrr::map_chr(const.SAMPLES_TEST, const.get_sample_uniform_chrAll)

const.SAMPLES_DNASEQ <- c(
  "F",
  "RE1",
  "RE2",
  "RE3"
)

const.get_orig_sample <- function(sample) {
  sample |>
  as.character() |>
  purrr::map_chr(
    function(x) {
      if (x %in% const.SAMPLES) {
        x
      } else if (x %in% const.SAMPLES_HOTSPOT) {
        const.SAMPLES[[which(x == const.SAMPLES_HOTSPOT)]]
      } else if (x %in% const.SAMPLES_SHUFFLE) {
        const.SAMPLES[[which(x == const.SAMPLES_SHUFFLE)]]
      } else if (x %in% const.SAMPLES_UNIFORM) {
        const.SAMPLES[[which(x == const.SAMPLES_UNIFORM)]]
      } else if (x %in% const.SAMPLES_UNIFORM_CHRALL) {
        const.SAMPLES[[which(x == const.SAMPLES_UNIFORM_CHRALL)]]
      } else if (x %in% const.SAMPLES_TEST) {
        const.SAMPLES[[which(x == const.SAMPLES_TEST)]]
      } else if (x %in% const.SAMPLES_TEST_SHUFFLE) {
        const.SAMPLES[[which(x == const.SAMPLES_TEST_SHUFFLE)]]
      } else if (x %in% const.SAMPLE_TEST_UNIFORM) {
        const.SAMPLES[[which(x == const.SAMPLE_TEST_UNIFORM)]]
      } else if (x %in% const.SAMPLE_TEST_UNIFORM_CHRALL) {
        const.SAMPLES[[which(x == const.SAMPLE_TEST_UNIFORM_CHRALL)]]
      } else if (x %in% const.SAMPLES_DNASEQ) {
        x
      } else if (x %in% const.SAMPLES_DNASEQ_DRAW) {
        const.SAMPLES[[which(x == const.SAMPLES_DNASEQ_DRAW)]]
      } else {
        stop(paste0("Invalid sample: ", x))
      }
    }
  )
}

const.ENZYME_F <- "F"
const.ENZYME_F_RE1 <- "F_RE1"
const.ENZYME_RE1 <- "RE1"
const.ENZYME_RE2 <- "RE2"
const.ENZYME_RE3 <- "RE3"
const.ENZYMES <- c(
  const.ENZYME_F,
  const.ENZYME_F_RE1,
  const.ENZYME_RE1,
  const.ENZYME_RE2,
  const.ENZYME_RE3
)

const.F_ENZYMES <- c(
  const.ENZYME_F,
  const.ENZYME_F_RE1
)

const.RE_ENZYMES <- c(
  const.ENZYME_RE1,
  const.ENZYME_RE2,
  const.ENZYME_RE3
)

const.CELLS <- c(
  "CD4T",
  "HCT116-WT",
  "HEK283T-WT",
  "HEK293T-RNASEH2A-KO-T3-8",
  "HEK293T-RNASEH2A-KO-T3-17",
  "hESC-H9",
  "WB-GTP"
)

const.CELL_CLASSES <- c("WT", "KO")

const.SAMPLE_NORMAL <- "normal"
const.SAMPLE_NORMAL_WT <- "normal_wt"
const.SAMPLE_NORMAL_KO <- "normal_ko"
const.SAMPLE_HOTSPOT <- "hotspot"
const.SAMPLE_SHUFFLE <- "shuffle"
const.SAMPLE_UNIFORM <- "uniform"
const.SAMPLE_UNIFORM_CHRALL <- "uniform_chrAll"
const.SAMPLE_TEST <- "test"
const.SAMPLE_TEST_SHUFFLE <- "test_shuffle"
const.SAMPLE_TEST_UNIFORM <- "test_uniform"
const.SAMPLE_TEST_UNIFORM_CHRALL <- "test_uniform_chrAll"
const.SAMPLE_NORMAL_ENZYME_F <- "normal_enzyme_f"
const.SAMPLE_UNIFORM_ENZYME_F <- "uniform_enzyme_f"
const.SAMPLE_UNIFORM_CHRALL_ENZYME_F <- "uniform_chrAll_enzyme_f"
const.SAMPLE_SHUFFLE_ENZYME_F <- "shuffle_enzyme_f"
const.SAMPLE_TEST_ENZYME_F <- "test_enzyme_f"
const.SAMPLE_DNASEQ <- "dnaseq"
const.SAMPLE_DNASEQ_DRAW <- "dnaseq_draw"
const.SAMPLE_DNASEQ_DRAW_ENZYME_F <- "dnaseq_draw_enzyme_f"

const.SAMPLE_TYPES <- c(
  const.SAMPLE_NORMAL,
  const.SAMPLE_NORMAL_WT,
  const.SAMPLE_NORMAL_KO,
  const.SAMPLE_HOTSPOT,
  const.SAMPLE_SHUFFLE,
  const.SAMPLE_UNIFORM,
  const.SAMPLE_UNIFORM_CHRALL,
  const.SAMPLE_TEST,
  const.SAMPLE_TEST_SHUFFLE,
  const.SAMPLE_TEST_UNIFORM,
  const.SAMPLE_TEST_UNIFORM_CHRALL,
  const.SAMPLE_NORMAL_ENZYME_F,
  const.SAMPLE_SHUFFLE_ENZYME_F,
  const.SAMPLE_UNIFORM_ENZYME_F,
  const.SAMPLE_UNIFORM_CHRALL_ENZYME_F,
  const.SAMPLE_TEST_ENZYME_F,
  const.SAMPLE_DNASEQ,
  const.SAMPLE_DNASEQ_DRAW,
  const.SAMPLE_DNASEQ_DRAW_ENZYME_F
)

const.SAMPLE_LIST <- list()
const.SAMPLE_LIST[[const.SAMPLE_NORMAL]] <- const.SAMPLES
const.SAMPLE_LIST[[const.SAMPLE_NORMAL_WT]] <- const.SAMPLES_WT
const.SAMPLE_LIST[[const.SAMPLE_NORMAL_KO]] <- const.SAMPLES_KO
const.SAMPLE_LIST[[const.SAMPLE_HOTSPOT]] <- const.SAMPLES_HOTSPOT
const.SAMPLE_LIST[[const.SAMPLE_SHUFFLE]] <- const.SAMPLES_SHUFFLE
const.SAMPLE_LIST[[const.SAMPLE_UNIFORM]] <- const.SAMPLES_UNIFORM
const.SAMPLE_LIST[[const.SAMPLE_UNIFORM_CHRALL]] <- const.SAMPLES_UNIFORM_CHRALL
const.SAMPLE_LIST[[const.SAMPLE_TEST]] <- const.SAMPLES_TEST
const.SAMPLE_LIST[[const.SAMPLE_TEST_SHUFFLE]] <- const.SAMPLES_TEST_SHUFFLE
const.SAMPLE_LIST[[const.SAMPLE_TEST_UNIFORM]] <- const.SAMPLES_TEST_UNIFORM
const.SAMPLE_LIST[[const.SAMPLE_TEST_UNIFORM_CHRALL]] <- const.SAMPLES_TEST_UNIFORM_CHRALL
const.SAMPLE_LIST[[const.SAMPLE_NORMAL_ENZYME_F]] <- const.SAMPLES_ENZYME_F
const.SAMPLE_LIST[[const.SAMPLE_SHUFFLE_ENZYME_F]] <- const.SAMPLES_SHUFFLE_ENZYME_F
const.SAMPLE_LIST[[const.SAMPLE_UNIFORM_ENZYME_F]] <- const.SAMPLES_UNIFORM_ENZYME_F
const.SAMPLE_LIST[[const.SAMPLE_UNIFORM_CHRALL_ENZYME_F]] <- const.SAMPLES_UNIFORM_CHRALL_ENZYME_F
const.SAMPLE_LIST[[const.SAMPLE_TEST_ENZYME_F]] <- const.SAMPLES_TEST_ENZYME_F
const.SAMPLE_LIST[[const.SAMPLE_DNASEQ]] <- const.SAMPLES_DNASEQ
const.SAMPLE_LIST[[const.SAMPLE_DNASEQ_DRAW]] <- const.SAMPLES_DNASEQ_DRAW
const.SAMPLE_LIST[[const.SAMPLE_DNASEQ_DRAW_ENZYME_F]] <- const.SAMPLES_DNASEQ_DRAW_ENZYME_F

const.get_sample_list <- function(sample_type) {
  if (sample_type %in% names(const.SAMPLE_LIST)) {
    const.SAMPLE_LIST[[sample_type]]
  } else {
    stop(paste0("Invalid sample type: ", sample_type))
  }
}

const.get_sample_info <- function(sample) {
  if(
    (sample %in% const.SAMPLES) ||
    (sample %in% const.SAMPLES_HOTSPOT) ||
    (sample %in% const.SAMPLES_SHUFFLE) ||
    (sample %in% const.SAMPLES_UNIFORM) ||
    (sample %in% const.SAMPLES_UNIFORM_CHRALL) ||
    (sample %in% const.SAMPLES_TEST) ||
    (sample %in% const.SAMPLES_TEST_SHUFFLE) ||
    (sample %in% const.SAMPLES_TEST_UNIFORM) ||
    (sample %in% const.SAMPLES_TEST_UNIFORM_CHRALL) ||
    (sample %in% const.SAMPLES_DNASEQ_DRAW)
  ) {
    const.SAMPLE_INFO[[const.get_orig_sample(sample)]]
  } else if (
    sample %in% const.SAMPLES_DNASEQ
  ) {
    const.SAMPLE_INFO_DNASEQ[[sample]]
  } else {
    stop("Invalid sample: ", sample)
  }
}

const.get_cell <- function(sample) {
  purrr::map_chr(
    sample,
    function(x) const.get_sample_info(x)[["cell"]]
  )
}

const.get_enzyme <- function(sample) {
  purrr::map_chr(
    sample,
    function(x) const.get_sample_info(x)[["enzyme"]]
  )
}

const.get_cell_class <- function(sample) {
  purrr::map_chr(
    sample,
    function(x) const.get_sample_info(x)[["cell_class"]]
  )
}

const.get_label <- function(sample) {
  sample |>
  purrr::map_chr(
    function(x) {
      if (
        (x %in% const.SAMPLES) ||
        (x %in% const.SAMPLES_HOTSPOT) ||
        (x %in% const.SAMPLES_SHUFFLE) ||
        (x %in% const.SAMPLES_UNIFORM) ||
        (x %in% const.SAMPLES_UNIFORM_CHRALL) ||
        (x %in% const.SAMPLES_TEST) ||
        (x %in% const.SAMPLES_TEST_SHUFFLE) ||
        (x %in% const.SAMPLES_TEST_UNIFORM) ||
        (x %in% const.SAMPLES_TEST_UNIFORM_CHRALL) ||
        (x %in% const.SAMPLES_DNASEQ_DRAW)
      ) {
        paste0(
          const.get_orig_sample(x),
          " ",
          const.get_cell(x),
          " ",
          const.get_enzyme(x),
          dplyr::case_when(
            (x %in% const.SAMPLES_HOTSPOT) ~ " (hotspot)",
            (x %in% const.SAMPLES_SHUFFLE) ~ " (shuffle)",
            (x %in% const.SAMPLES_UNIFORM) ~ " (uniform)",
            (x %in% const.SAMPLES_UNIFORM_CHRALL) ~ " (uniform chrAll)",
            (x %in% const.SAMPLES_TEST) ~ " (test)",
            (x %in% const.SAMPLES_TEST_SHUFFLE) ~ " (test shuffle)",
            (x %in% const.SAMPLES_TEST_UNIFORM) ~ " (test uniform)",
            (x %in% const.SAMPLES_TEST_UNIFORM_CHRALL) ~ " (test uniform chrAll)",
            (x %in% const.SAMPLES_DNASEQ_DRAW) ~ " (DNA-seq random)",
            TRUE ~ ""
          )
        )
      } else if (
        (x %in% const.SAMPLES_DNASEQ)
      ) {
        paste0(
          const.get_cell(x),
          " ",
          x,
          dplyr::case_when(
            x %in% const.SAMPLES_DNASEQ ~ " (DNA-seq)",
            TRUE ~ ""
          )
        )
      } else {
        stop("Invalid sample: ", x)
      }
    }
  )
}

const.get_sample_label_factor <- function(sample, sample_list) {
  sample |>
  factor(levels = sample_list) |>
  forcats::fct_relabel(const.get_label)
}

const.get_dnaseq_sample <- function(sample) {
  "F"
}

const.BED_COLS <- c("chrom", "start", "end", "name", "score", "strand")
