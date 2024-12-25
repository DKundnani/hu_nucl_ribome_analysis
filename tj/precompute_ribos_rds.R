# Convert BED files to RDS for easier usage with R

ribos_rds.do_main <- function(sample_list = const.SAMPLES, overwrite = FALSE) {
  purrr::walk(
    sample_list,
    function(sample) {
      if (
        const.CHROMS |>
        purrr::map_lgl(\(x) file.exists(fn.ribos_rds(sample = sample, chrom = x))) |>
        all()
      ) {
        utils.log("Already exists:", sample)
        return(invisible())
      }

      utils.read_tsv(fn.ribos_tsv(sample)) |>
      dplyr::filter(chrom %in% const.CHROMS) |>
      dplyr::arrange(chrom, strand, pos) |>
      x => {
        GenomicRanges::GRanges(
          seqnames = x[["chrom"]],
          ranges = IRanges::IRanges(
            start = x[["pos"]],
            end = x[["pos"]]
          ),
          strand = x[["strand"]],
          sample = S4Vectors::Rle(sample, nrow(x)),
          ribo_nuc = x[["ribo_nuc"]],
          count = x[["count"]]
        )
      } |>
      x => GenomicRanges::split(
        x,
        factor(GenomicRanges::seqnames(x), const.CHROMS) # splitting on factor always has all levels
      ) |>
      purrr::iwalk(
        function(x, chrom) {
          utils.write_rds(
            if (NROW(x) == 0) {
              utils.log("No data for", sample, chrom)
              NULL
            } else {
              x
            },
            fn.ribos_rds(sample = sample, chrom = chrom)
          )
        }
      )
    }
  )
}
