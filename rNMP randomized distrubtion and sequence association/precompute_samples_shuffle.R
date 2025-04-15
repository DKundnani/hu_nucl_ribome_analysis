# Shuffle the true samples to get shuffled samples
samples_shuffle.do_main <- function(
  sample_list = c(const.SAMPLES, const.SAMPLES_TEST),
  chrom_list = const.CHROMS,
  overwrite = FALSE
) {
  set.seed(123321) # For reproducibility

  gaps <- utils.get_hg38_gaps() # For removing gaps

  purrr::walk(
    sample_list,
    function(sample) {
      sample_shuffle <- const.get_sample_shuffle(sample)

      purrr::walk(
        chrom_list,
        function(chrom) {
          file_out <- fn.ribos_rds(sample = sample_shuffle, chrom = chrom)
          if (file.exists(file_out) && !overwrite) {
            utils.log("Already exists:", file_out)
            return(invisible())
          }

          granges <- utils.load_ribos_rds(sample = sample, chrom = chrom)
          chrom_size <- utils.get_chrom_sizes(chrom)

          start_shuffle <- sample(chrom_size, 2 * NROW(granges))
          granges_shuffle <- GenomicRanges::GRanges(
            seqnames = chrom,
            ranges = IRanges::IRanges(
              start = start_shuffle,
              end = start_shuffle
            )
          )
          granges_shuffle <- granges_shuffle[!(granges_shuffle %over% gaps),] # remove gaps
          if (NROW(granges_shuffle) < NROW(granges)) {
            # This should really never happen,
            # 2 * NROW(granges) should be enough even after removing gaps
            utils.log("Shuffled granges too small:", sample, chrom)
          }
          granges_shuffle <- granges_shuffle[seq_len(NROW(granges)),] # Trim to size of granges

          GenomicRanges::strand(granges_shuffle) <- strand(granges)
          granges_shuffle$sample <- S4Vectors::Rle(
            sample_shuffle,
            NROW(granges_shuffle)
          )
          granges_shuffle$count <- sample(granges$count)
          granges_shuffle <- granges_shuffle[order(granges_shuffle),]
          granges_shuffle$ribo_nuc <- (
            Biostrings::getSeq(
              BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
              granges_shuffle
            ) |>
            as.character()
          )
          # remove "N" in nucs, although it really should not be there
          # since we already removed gaps
          granges_shuffle <- granges_shuffle[granges_shuffle$ribo_nuc %in% const.RIBO_NUCS,]
          utils.write_rds(granges_shuffle, file_out)
        }
      )
    }
  )
}
