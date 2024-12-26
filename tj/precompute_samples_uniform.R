# Do uniform random sampling to generate new samples.
# This script has two sampling versions:
# 1. Sample each chromosome separately.
# 2. Sample all chromosomes at once.

# Sample the desired number of ribos uniformly from a single chromosome.
samples_uniform.sample_chrom <- function(sample_uniform, chrom, n_ribos) {
  gaps <- utils.get_hg38_gaps() # For removing gaps
  chrom_size <- utils.get_chrom_sizes(chrom)

  pos_uniform <- sample(chrom_size, 2 * n_ribos, replace = TRUE)
  granges_uniform <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(
      start = pos_uniform,
      end = pos_uniform
    )
  )
  granges_uniform <- granges_uniform[!(granges_uniform %over% gaps),] # remove gaps
  if (NROW(granges_uniform) < n_ribos) {
    # This should really never happen,
    # 2 * n_ribos should be enough even after removing gaps
    stop("Sampled granges too small:", sample, chrom)
  }
  granges_uniform <- granges_uniform[seq_len(n_ribos),] # Trim to size of granges
  granges_uniform$count <- 1 # 1 ribo at each position

  # combine the counts of the overlapping ribos
  count <- tapply(
    granges_uniform$count,
    GenomicRanges::start(granges_uniform),
    sum
  )
  pos <- as.integer(names(count))
  granges_uniform <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(
      start = pos,
      end = pos
    ),
    strand = sample(const.STRANDS, length(count), replace = TRUE), # random strands
    sample = S4Vectors::Rle(sample_uniform, length(count)),
    count = as.numeric(count)
  )
  granges_uniform <- granges_uniform[order(granges_uniform),]
  
  # Check if there are any Ns in the nucleotide sequence.
  # There should not be any since we already removed gaps.
  granges_uniform$ribo_nuc <- ( # get the nucleotide sequence
    Biostrings::getSeq(
      BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
      granges_uniform
    ) |>
    as.character()
  )
  if (any(granges_uniform$ribo_nuc %in% "N")) {
    stop("N in nucs:", sample, chrom)
  }

  granges_uniform
}

# Make random samples that have the same number of ribos as the
# original samples in each chromosome.
samples_uniform.do_main <- function(
  sample_list = c(const.SAMPLES, const.SAMPLES_TEST),
  chrom_list = const.CHROMS,
  overwrite = FALSE
) {
  set.seed(54765) # For reproducibility

  purrr::walk(
    sample_list,
    function(sample) {
      sample_uniform <- const.get_sample_uniform(sample)

      purrr::walk(
        chrom_list,
        function(chrom) {
          file_out <- fn.ribos_rds(sample = sample_uniform, chrom = chrom)
          if (file.exists(file_out) && !overwrite) {
            utils.log("Already exists:", file_out)
            return(invisible())
          }

          granges <- utils.load_ribos_rds(sample = sample, chrom = chrom)
          n_ribos <- sum(granges$count)
          utils.write_rds(
            samples_uniform.sample_chrom(sample_uniform, chrom, n_ribos),
            file_out
          )
        }
      )
    }
  )
}

# Make random samples that have the same number of ribos as the
# original samples in the whole genome but not necessarily in each chromosome.
samples_uniform.do_main_chrAll <- function(
  sample_list = c(const.SAMPLES, const.SAMPLES_TEST),
  overwrite = FALSE
) {
  set.seed(24534) # For reproducibility

  purrr::walk(
    sample_list,
    function(sample) {
      sample_uniform <- const.get_sample_uniform_chrAll(sample)

      file_exists <- purrr::map_lgl(
        const.CHROMS,
        \(chrom) file.exists(fn.ribos_rds(sample = sample_uniform, chrom = chrom))
      )
      if (all(file_exists) && !overwrite) {
        utils.log("Already exists:", sample_uniform)
        return(invisible())
      }
      
      granges <- utils.load_ribos_rds(sample = sample, chrom = const.CHROM_ALL)

      # Get the size of each chromosome and its gaps.
      gap_size <- (
        utils.get_hg38_gaps() |>
        utils.granges_to_tibble() |>
        dplyr::group_by(chrom) |>
        dplyr::summarize(gap_size = sum(width)) |>
        x => setNames(x[["gap_size"]], x[["chrom"]])
      )
      chrom_size <- utils.get_chrom_sizes(const.CHROMS)

      n_ribos <- sum(granges$count) # Number of ribos in the whole genome.

      # First stage: sample the chromosomes themselves based on their non-gap length.
      # We use non-gap length because we are not allowed to sample ribos from gaps.
      n_ribos_chrom <- ( # Number of ribos in each chromosome.
        sample(
          const.CHROMS,
          size = n_ribos,
          prob = chrom_size[const.CHROMS] - gap_size[const.CHROMS],
          replace = TRUE
        ) |>
        table()
      )

      # Second stage: sample from each chromosome/strand uniformly.
      purrr::walk2(
        const.CHROMS,
        n_ribos_chrom[const.CHROMS],
        function(chrom, n_ribos) {
          utils.write_rds(
            samples_uniform.sample_chrom(sample_uniform, chrom, n_ribos),
            fn.ribos_rds(sample = sample_uniform, chrom = chrom)
          )
        }
      )
    }
  )
}
