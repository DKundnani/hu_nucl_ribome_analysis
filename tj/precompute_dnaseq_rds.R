# Converts DNA-seq BED files to RDS for easier usage with R

dnaseq_rds.do_main <- function(
  sample_list = const.SAMPLES_DNASEQ,
  chrom_list = const.CHROMS,
  overwrite = FALSE
) {
  purrr::pwalk(
    tidyr::expand_grid(
      sample = sample_list,
      chrom = chrom_list
    ),
    function(sample, chrom) {
      file_out <- fn.dnaseq_rds(
        sample = sample,
        assembly = "hg38",
        chrom = chrom
      )
      if (!overwrite && file.exists(file_out)) {
        utils.log("Already exists:", file_out)
        return(invisible())
      }

      fn.dnaseq_bed(
        enzyme = sample,
        assembly = "hg38",
        chrom = chrom
      ) |>
      utils.read_tsv(
        col_names = c(
          "chrom",
          "start",
          "end",
          "covered",
          "count",
          "strand"
        ),
        col_types = "cnnlnc"
      ) |>
      dplyr::select(chrom, strand, start, end, count) |>
      utils.format_bed_as_table() |>
      dplyr::filter(chrom %in% const.CHROMS) |>
      dplyr::arrange(chrom, strand, start) |>
      x => {
        GenomicRanges::GRanges(
          seqnames = factor(chrom, const.CHROMS),
          ranges = IRanges::IRanges(
            start = x[["start"]],
            end = x[["end"]]
          ),
          strand = x[["strand"]],
          sample = S4Vectors::Rle(sample, nrow(x)),
          count = x[["count"]]
        )
      } |>
      x => (if (NROW(x) == 0) NULL else x) |>
      utils.write_rds(file_out)
    }
  )
}
