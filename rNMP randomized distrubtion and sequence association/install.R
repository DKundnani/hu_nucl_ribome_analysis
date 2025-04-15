## Copied/modified from ../ballard/0-install.R 

## Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

## Install packages
BiocManager::install(
  c(
    "rtracklayer", # importing/exporting bed files (and others)
    "BSgenome.Hsapiens.UCSC.hg38", # reference human genome
    "BSgenome.Hsapiens.UCSC.hg38.masked", # reference human genome with masks
    "GenomicRanges", # GRanges objects
    "IRanges", # IRanges objects
    "plyranges", # for manipulating GRanges objects
    "Biostrings", # utilies for DNA/RNA sequences (e.g., couting kmers)
    "annotatr" # for getting annotations such as genes, CpG islands, etc.
  )
)

## General utilities
install.packages("tidyverse")
install.packages("xtable")
install.packages("lubridate")
install.packages("gtools")
install.packages("pheatmap")

# For "wasserstein1d" function
install.packages("transport")

# For dynamic time warp distance
install.packages("dtw")
