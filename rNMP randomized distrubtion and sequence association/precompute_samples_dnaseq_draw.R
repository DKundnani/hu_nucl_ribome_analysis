# Use the DNA-seq samples as a background distribution to draw random samples
samples_dnaseq_draw.do_main <- function(
  sample_list = const.SAMPLES,
  chrom_list = const.CHROMS,
  overwrite = FALSE
) {
  set.seed(123321)

  purrr::walk(
    chrom_list,
    function(chrom) {
      data_dnaseq <- (
        utils.load_dnaseq_rds(
          sample = "F",
          chrom = chrom,
          assembly = "hg38"
        ) |>
        utils.granges_to_tibble() |>
        dplyr::filter(count > 0) |>
        dplyr::nest_by(strand)
      )

      purrr::walk(
        sample_list,
        function(sample) {
          sample_dnaseq_draw <- const.get_sample_dnaseq_draw(sample)
          file_out <- fn.ribos_rds(
            sample = sample_dnaseq_draw,
            chrom = chrom
          )
          if (!overwrite && file.exists(file_out)) {
            utils.log("Already exists:", file_out)
            return(invisible())
          }

          ribo_data <- (
            fn.ribos_rds(
              sample = sample,
              chrom = chrom
            ) |>
            utils.read_rds()
          )

          if (is.null(ribo_data)) {
            utils.log("No ribos for", sample, chrom)
            utils.write_rds(NULL, file_out)
            return(invisible())
          }

          data_count_total <- (
            ribo_data |>
            utils.granges_to_tibble() |>
            dplyr::group_by(strand) |>
            dplyr::summarize(count_total = sum(count))
          )

          dplyr::inner_join(
            data_dnaseq,
            data_count_total,
            by = "strand"
          ) |>
          purrr::pmap(
            function(strand, data_dnaseq, count_total) {
              utils.log(sample, chrom, strand)
              widths <- data_dnaseq[["end"]] - data_dnaseq[["start"]] + 1
              weights <- data_dnaseq[["count"]] * widths

              # Two stage sampling:
              # 1. idx is the indices of the ranges we sample from
              # 2. offsets are the offsets within those ranges
              idx <- sample.int(length(weights), count_total, replace = TRUE, prob = weights)
              offsets <- widths[idx] |> purrr::map_int(sample.int, 1) # for each width, w, select an integer in [1, w]
              
              start <- data_dnaseq[["start"]][idx] + offsets - 1 # the final positions
              count <- table(start) # count the number of times each position was selected
              start <- as.integer(names(count)) # the final unique positions
              count <- as.integer(count) # the final counts
              end <- start
              ranges <- GenomicRanges::GRanges(
                seqnames = rep(chrom, length(start)),
                ranges = IRanges::IRanges(start, end),
                strand = rep(strand, length(start)),
                sample = S4Vectors::Rle(sample_dnaseq_draw, length(start)),
                count = count
              )
              ribo_nuc <- (
                Biostrings::getSeq(
                  BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                  ranges
                ) |>
                as.character()
              )
              ranges$ribo_nuc <- ribo_nuc
              ranges$count
              ranges <- ranges[ribo_nuc %in% const.RIBO_NUCS]
              ranges <- sort(ranges)
              ranges
            }
          ) |>
          x => do.call(plyranges::bind_ranges, x) |>
          x => {
            if (NROW(x) == 0) {
              NULL
            } else {
              x
            }
          } |>
          utils.write_rds(file_out)
        }
      )
    }
  )
}
