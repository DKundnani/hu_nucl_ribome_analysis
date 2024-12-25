init.init <- function(
  no_repeat = FALSE, # Controls whether we use the repeat-removed data or not.
  no_title = FALSE, # Controls whether all plots have the title removed or not.
  plot_file_type = "png", # "png" or "pdf".
  plot_dpi = 300, # dots per inch.
  plot_sample_label_margin = 6, # Padding to use in plots when the sample name is shown since some are very long; in inches.
  plot_title_margin = 6, # Padding to use in plots when the title is shown since some are very long; in inches.
  base_font_size = 18, # Base font size to use in plots.
  output_extra_suffix = "", # In case we want to isolate some analyses, use a suffix to distinguish them.
  plot_extra_suffix = "", # In case we want to isolate some plots, use a suffix to distinguish them.
  data_dir_base = Sys.getenv("DATA"), # Base directory for data; checks the environment variable "DATA" first.
  log_enable = TRUE, # Whether to print log message or not.
  test = FALSE # Whether to run analyses in test mode or not.
) {
  # Set up the settings global variables.
  assign("settings.NO_REPEAT", no_repeat, envir = .GlobalEnv)
  assign("settings.NO_TITLE", no_title, envir = .GlobalEnv)

  assign("settings.PLOT_FILE_TYPE", plot_file_type, envir = .GlobalEnv)
  assign("settings.PLOT_DPI", plot_dpi, envir = .GlobalEnv)
  assign("settings.PLOT_SAMPLE_LABEL_MARGIN", plot_sample_label_margin, envir = .GlobalEnv)
  assign("settings.PLOT_TITLE_MARGIN", plot_title_margin, envir = .GlobalEnv)

  assign("settings.BASE_FONT_SIZE", base_font_size, envir = .GlobalEnv)

  assign("settings.OUTPUT_EXTRA_SUFFIX", output_extra_suffix, envir = .GlobalEnv)
  assign("settings.PLOT_EXTRA_SUFFIX", plot_extra_suffix, envir = .GlobalEnv)

  assign("settings.DATA_DIR_BASE", data_dir_base, envir = .GlobalEnv)
  if (settings.DATA_DIR_BASE == "") {
    assign("settings.DATA_DIR_BASE", ".", envir = .GlobalEnv) # default to current directory
  }

  assign("settings.LOG_ENABLE", log_enable, envir = .GlobalEnv)

  assign("settings.TEST", test, envir = .GlobalEnv)

  # Load external libraries
  library(tidyverse)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  library(GenomicRanges)
  library(plyranges)
  library(IRanges)
  library(Biostrings)
  library(xtable)
  library(lubridate)
  library(gtools)
  library(rtracklayer)

  # Load utilities
  source("channagiri/util_constants.R")
  source("channagiri/util_file_names.R")
  source("channagiri/util_kmer.R")
  source("channagiri/util_functions.R")
}