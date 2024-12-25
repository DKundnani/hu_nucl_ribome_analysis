samples_test.do_main <- function(
  sample_list = const.SAMPLES,
  chrom_list = const.CHROMS,
  test_size = 1000, # Number of rows to sample from each strand
  overwrite = FALSE
) {
  set.seed(123) # For reproducibility
  purrr::pwalk(
    tidyr::expand_grid(
      sample = sample_list,
      chrom = chrom_list
    ),
    function(sample, chrom) {
      sample_test <- const.get_sample_test(sample)
      file_out <- fn.ribos_rds(sample = sample_test, chrom = chrom)
      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }
      granges <- utils.load_ribos_rds(sample = sample, chrom = chrom)
      granges_plus <- granges[GenomicRanges::strand(granges) == "+",]
      granges_minus <- granges[GenomicRanges::strand(granges) == "-",]
      granges_test <- c(
        granges_plus[sample(NROW(granges_plus), min(NROW(granges_plus), test_size)),],
        granges_minus[sample(NROW(granges_minus), min(NROW(granges_minus), test_size)),]
      )
      granges_test <- granges_test[order(granges_test),]
      utils.write_rds(granges_test, file_out)
    }
  )
}
