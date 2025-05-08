# rNMP randomzed distrubtion and sequence association

## Installation

R version 4.2.2 or greater must be installed for these scripts to work. Please run the `install.R` script. This script will install the necessary prerequisites for running the other analyses. Installation may take a long time and a large amount of hard drive space because several [Bioconductor](https://www.bioconductor.org/) packages, including the entire [hg38](https://genome.ucsc.edu/cgi-bin/hgGateway) (GrCh38) reference genome, are downloaded. The version of Bioconductor used while developing these scripts is 3.16, though other versions may also work.

## Input data

* The BED files of mapped ribos locations must be saved to the directory `data/ribo`. This directory must contain the files: `FS185.tsv`, `FS186.tsv`, `FS187.tsv`, `FS188.tsv`, `FS189.tsv`, `FS193.tsv`, `FS195.tsv`, `FS196.tsv`, `FS197.tsv`, `FS198.tsv`, `FS199.tsv`, `FS201.tsv`, `FS203.tsv`, `FS300.tsv`, `FS301.tsv`, `FS303.tsv`, `FS305.tsv`, `FS306.tsv`, `FS307.tsv`, `FS309.tsv`, `FS310.tsv`, `FS326.tsv`, `FS327.tsv`, `FS329.tsv`, `FS331.tsv`, `FS333.tsv`, `FS391.tsv`, `FS392.tsv`, and `FS393.tsv`. Each `.tsv` file must be a tab-separated table with columns `chrom`, `strand`, `pos`, `count`, and `ribo_nuc`. The data is currently at
  [rNMP_Project_data_2025_04_21/channagiri/data/ribos](rNMP_Project_data_2025_04_21/channagiri/data/ribos)
  <!--[https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/rNMP_libraries_2024_01/full_tsv/](https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/rNMP_libraries_2024_01/full_tsv/).--!>
* The DNA-seq data of background read coverage on the reference genome must be saved to the `data/dnaseq` directory. Because the DNA-seq data is large, the data must already be split into separate chromosomes. The following subdirectories must be present: `data/dnaseq/cd4t_dsF_hg38_coverage`, `data/dnaseq/cd4t_RE1_hg38_coverage`, `data/dnaseq/cd4t_RE2_hg38_coverage`, `data/dnaseq/cd4t_RE3_hg38_coverage`. Each subdirectory must have (at least) the 23 files: `chr1.bed`, `chr2.bed`, `chr3.bed`, `chr4.bed`, `chr5.bed`, `chr6.bed`, `chr7.bed`, `chr8.bed`, `chr9.bed`, `chr10.bed`, `chr11.bed`, `chr12.bed`, `chr13.bed`, `chr14.bed`, `chr15.bed`, `chr16.bed`, `chr17.bed`, `chr18.bed`, `chr19.bed`, `chr20.bed`, `chr21.bed`, `chr22.bed`, `chrX.bed`. The data is currently at [rNMP_Project_data_2025_04_21/DNAseq_coverage/full_split/](rNMP_Project_data_2025_04_21/DNAseq_coverage/full_split/)
  <!--[https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/DNAseq_coverage/full_split/](https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/DNAseq_coverage/full_split/).--!>
* The hg38 gaps (i.e., regions filled with Ns) BED file must be saved to `data_common/hg38_gaps.bed`. The file is currently at
[rNMP_Project_data_2025_04_21/UCSC_hg38/hg38_gaps.bed](rNMP_Project_data_2025_04_21/UCSC_hg38/hg38_gaps.bed)
<!--[https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/UCSC_hg38/hg38_gaps.bed](https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/UCSC_hg38/hg38_gaps.bed).--!>


An easier way to download the data files is to mirror the subdirectories `data` and `data_common` at [https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/rNMP_Hotspots/channagiri](https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/rNMP_Hotspots/channagiri). Another way is to download the compressed archive with all the files at [https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/rNMP_Hotspots_data](https://knot.math.usf.edu/safe/users/channagiri/rNMP_Project/rNMP_Hotspots_data) (if there are multiple archives, please use the latest).

## Settings

Please see `settings.R` for important settings related to the analysis and output.

## Utility scripts

No action is needed here.

* `util_constants.R`: Useful constants decribing the samples, chromosome, etc.
* `util_file_names.R`: Functions for forming output file names.
* `util_functions.R`: General utility functions that are used in different analyses.
* `util_kmer.R`: Utiliy functions/constants for working with kmers.

## Initialization

The `init.R` script loads the necessary libraries and files for running the precomputations and analyses. No action is needed here.

## Precomputation

Once the BED files are saved, please run the `run_precompute.R` script to do the precomputation. This will do the following.

* `weights_binned.do_main()`: Get the number of A, C, G, T nucleotides (that is, non-N nucleotides) in each window of the reference genome for different window sizes.
* `kmer_binned.do_main()`: Get the number of each kmer (length 1, 2, or 3) in each window of the reference genome, for different window sizes. 
* `ribos_rds.do_main()`: Split the [BED](https://useast.ensembl.org/info/website/upload/bed.html) files by chromosome, convert them to [GRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) objects, and write them to [RDS format](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html), for quicker loading.
* `ribos_binned.do_main()`: Compute the number of ribos in each window of the reference genome, for different window sizes.
* `samples_test.do_main()`: Create test samples, which are just small subsets of the full data sets. These are used for testing analyses.
* `samples_shuffle.do_main()`: Create shuffle samples, which may be used as control data sets.
* `samples_uniform.do_main()`: Create samples by drawing "ribos" uniformly at random within each chromosome. These may be used as control data sets.
* `samples_uniform_chrAll.do_main()`: Create samples by drawing "ribos" uniformly at random within the entire genome. These may be used as control data sets.
* `dnaseq_rds.do_main()`: Convert the DNA-seq [BED](https://useast.ensembl.org/info/website/upload/bed.html) files to [GRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) objects, and write them to [RDS format](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html), for quicker loading.
* `dnaseq_binned.do_main()`: Compute the DNA-seq count in each window of the reference genome, for different window sizes. The count is in read*base units, which is the sum of the number of reads overlapping each base for all bases in the window. Dividing this count by the number of non-N nucleotides in the window will give you the mean coverage of that window.
* `samples_dnaseq_draw.do_main()`: Create DNA-seq draw samples, which may be used as control data sets. These samples are generated by using the DNA-seq distribution as the background distribution.
* `kmer_heatmaps.do_precompute()`: Precompute the kmer heatmap data, which counts the normalized frequency of kmers within a 50 bp window of each ribo in the sample (multiple ribos in a single location are counted multiple times).

The individual precomputations steps are implemented in the R scripts beginning with `precompute_`.

## Analysis

Once the precomputations is finished, please run the analysis with the `run_analysis.R` script. This will do the following:

* `summary.do_main()`: Create summary tables describing the distribution of rNMPs in each sample and chromosome.
* `kmer_heatmaps.do_main()`: Create heatmap figures showing the frequency (or some normalization thereof) of kmers within a certain window of the rNMPs.
* `kmer_correlation.do_main()`: Create heatmap figures showing the correlation of ribo frequencies (or some normalization thereof) with kmer frequencies. This uses binned data.
* `kmer_hist_offset.do_main()`: Create line graph figures showing the frequency (or some normalization thereof) at different offsets from the ribo locations. This is resource intensive and may take a lot of time to complete.
* `distribution.do_main()`: Create figures comparing the distribution of rNMPs in the regular samples with the control samples (e.g., DNA-seq draw or shuffle samples).
* `line_graph.do_main()`: Create line graph figures showing the frequency (or some normalization thereof) of rNMPs along each chromosome. Also creates heatmaps showing the pairwise correlations (or other pairwise quantity) of samples. This uses binned data.
* `neighborhood_counts.do_main()`: Creates figures showing the distribution of "neighborhood counts", which are the the number of rNMPs within a neighborhood of each ribo. This is intended to allow judging how much clustering there is at rNMPs locations.
* `between_ribo_dist.do_main()`: Creates figures showing the distribution of "between ribo distances", which is the distance between consecutive rNMPs. This is intended to allow judging how much clustering there is at rNMPs locations.
* `dnaseq_analysis.do_main()`: Creates figures analyzing the DNA-seq data. This includes the distribution of the DNA-seq data, their correlation with ribo samples, etc. This mostly uses binned data.
* `overlap_clusters.do_main()`: Creates figures showing how many "overlap clusters" of ribo locations there are with a certain window width. An overlap cluster is a connected component obtained when using the relation "within D base pairs of each other" on the ribo locations, for some D.

To simply test the analyses for development purposes, set `settings.TEST = TRUE` in `settings.R`. The individual analyses are implemented in the R scripts beginning with `analysis_`.

## Output

The directories with the following prefixes contain the output.

* `data_output`: Contains intermediate processing tables and analysis output.
* `plot`: Contains output figures.

The suffixes may indicate different settings that have been set in `settings.R` (e.g., `plot_png` will have PNG files, `plot_pdf` will have PDF files, etc.,).

## Latex

The following directories contain LaTeX source and output files explaining some of the analyses.

* `latex`: LaTeX source files. When building the LaTeX files, plots must sometimes be generated first by running the appropriate R files in the directory.
* `latex_pdf`: PDF files generated from the LaTeX sources. Not all the PDFs are complete or even relevant anymore. Please always refer to the code for the most current information.

## File naming conventions

Output figure and table files are named according to the parameters used in their creation. These may include, for example, the window size used for binning, the sample name, the chromosome, etc. The following in an example file name.

```
matrix_uniform_100kb_corr_pearson_chr1_-_nucAll_divideNonGap.png
```

* `matrix` indicates that this figure is a pairwise correlation matrix.
* `uniform` indicates that this figure used the "uniform" control sample (see below for more details).
* `100kb` indicates that this figure used a window size of 100,000 bp.
* `corr_pearson` indicates that this figure used the Pearson correlation coefficient.
* `chr1` indicates that this figure is for chromosome 1.
* `-` indicates that this figure is for the minus strand.
* `nucAll` indicates that this figure used data for all ribo nucleotides (rA, rC, rG, rU).
* `divideNonGap` indicates that this figure is normalized by dividing by the number of proper (aka "non-gap") nucleotides (A, C, G, T). A "gap" nucleotide by constast is the letter N where no rNMP data may exist.

## Control samples

In the file names, the type of control or experimental sample are indicated by the following:

* `normal`: The experimental rNMP data (not a control).
* `uniform`: Control sample that is a uniform distribution of rNMPs across the genome.
* `shuffle`: Control sample that is a shuffle of the original sample.
* `dnaseq_draw`: Control sample that is a random sampes using the Fragmentase DNA-seq distribution as a background.

## Normalization methods

In the file names, the normalization method is indicated by the following.

* `divideNonGap`: The normalization is done by dividing the number of proper (aka "non-gap") nucleotides (A, C, G, T). A "gap" nucleotide by constrast is the letter N where no rNMP data may exist.
* `width`: The normalization is done by dividing by the width of the window.
* `none`: No normalization is done.

Please see the `utils.load_ribos_binned()` method in the `util_functions.R` script for more details.

## Window sizes

In the file names for binned analyses, the window size is indicated by the following:

* `1mb`: The window size is 1,000,000 bp.
* `100kb`: The window size is 100,000 bp.
* `10kb`: The window size is 10,000 bp.
* `1kb`: The window size is 1,000 bp.

## rNMP nucleotides

In the file names, the rNMP ribonucleotides are indicated by the following.

* `nucAll`: Data from all ribo nucleotides (rA, rC, rG, rU) is used.
* `A`: Data from only rA is used.
* `C`: Data from only rC is used.
* `G`: Data from only rG is used.
* `T`: Data from only rU is used.

## Strand

In the file names, the strand is indicated by the following.

* `+-`: Data from both the plus and minus strands are used.
* `+`: Data from the plus strand is used.
* `-`: Data from the minus strand is used.

## Chromosome

In the file names, the chromosome is indicated by the following.

* `chrAll`: Data from all chromosomes is used.
* `chr1`: Data from Chromosome 1 is used.
* `chr2`: Data from Chromosome 2 is used.
* `chr3`: Data from Chromosome 3 is used.
* `chr4`: Data from Chromosome 4 is used.
* `chr5`: Data from Chromosome 5 is used.
* `chr6`: Data from Chromosome 6 is used.
* `chr7`: Data from Chromosome 7 is used.
* `chr8`: Data from Chromosome 8 is used.
* `chr9`: Data from Chromosome 9 is used.
* `chr10`: Data from Chromosome 10 is used.
* `chr11`: Data from Chromosome 11 is used.
* `chr12`: Data from Chromosome 12 is used.
* `chr13`: Data from Chromosome 13 is used.
* `chr14`: Data from Chromosome 14 is used.
* `chr15`: Data from Chromosome 15 is used.
* `chr16`: Data from Chromosome 16 is used.
* `chr17`: Data from Chromosome 17 is used.
* `chr18`: Data from Chromosome 18 is used.
* `chr19`: Data from Chromosome 19 is used.
* `chr20`: Data from Chromosome 20 is used.
* `chr21`: Data from Chromosome 21 is used.
* `chr22`: Data from Chromosome 22 is used.
* `chrX`: Data from Chromosome X is used.

## Between ribo distance analysis names

The following are parts of file names in the between ribo distance analysis.

* `cdf`: The cumulative distribution function.
* `pdf`: The probability density function.

## Distribution analysis names

The following are parts of file names in the distribution analysis.

* The first part of the figure names shows the pair of samples being compared. For example `normal_uniform` means the experimental rNMP sample is compared with the uniform control sample.
* `density`: Histogram of the distribution of ribo counts in windows comparing the two samples types.
* `pp`: Probability-probability plot of the distribution of ribo counts in windows comparing the two sample types.
* `qq`: Quantile-quantile plot of the distribution of ribo counts in windows.

## DNA-seq analysis names

The following are parts of file names in the DNA-seq analysis.

* `matrix`: Heatmap of the pairwise correlation of ribo samples with DNA-seq samples. The type of correlation is written into the file name (`pearson`, `spearman`, or `kendall`).
* `line_graph`: Line graph of the frequency of DNA-seq along the chromosome.
* `line_graph_ribo`: Line graph of the frequency of DNA-seq vs ribos along the chromosome.
* `noOutlier`: The DNA-seq samples are filtered to remove outliers. An outlier is a window whose count is more than 4 median-absolute-deviations (MADs) above the median count. Please see the `utils.load_dnaseq_binned()` function in the `util_functions.R` script for more details.
* `outlier`: The DNA-seq samples are not filtered to remove outliers.

## Kmer correlation analysis names

The following are parts of file names in the kmer correlation analysis.

* `matrix`: Heatmap of the pairwise correlation of rNMP samples with kmer samples. The type of correlation is written into the file name (`pearson`, `spearman`, or `kendall`).
* `line`: Line graph showing the frequency of rNMPs vs background kmer frequencies along the chromosome.
* `scatter`: Scatter plot showing the frequency of rNMPs vs background kmer frequencies along the chromosome.
* `none`, `ratio`, `sum1`: The normalization method used for the kmer frequencies. `none` means no normalization is done. `ratio` means the kmer frequencies are divided by the background kmer frequency. `sum1` means the kmer frequencies are divided by the background kmer frequencies and normalized again to sum to one within each sample.
* `rank`, `sum`, `mean`, `none`: The second normalization scheme used after the first (above). `rank` converts the frequencies to ranks. `sum` divides the frequencies by the sum of the frequencies. `mean` divides the frequencies by the mean of the frequencies. This is done to allow easier comparison of the kmer frequencies with the rNMP frequencies.
* The size of the kmers used is written into the file name.
* `all` indicates that the figure is generated with all the kmers of the given size. If a subset is used, the subset is written into the file name instead.

## Kmer heatmap analysis names

The following are parts of file names in the kmer heatmap analysis.

* `up`, `down`, `both`: Whether the figure uses data from the upstream window, downstream window, or both windows around the rNMP location.
* The radius of the window around the rNMP location used to compute the kmer frequencies is written into the file name.
* `none`, `ratio`, `sum1`: The normalization method used for the kmer frequencies. `none` means no normalization is done. `ratio` means the kmer frequencies are divided by the background kmer frequency. `sum1` means the kmer frequencies are divided by the background kmer frequencies and normalized again to sum to one within each sample.
* The size of the kmer is written into the file name.

## Kmer offset histogram analysis names

The following are parts of file names in the kmer offset histogram analysis.

* `none`, `ratio`, `sum1`: The normalization method used for the kmer frequencies. `none` means no normalization is done. `ratio` means the kmer frequencies are divided by the background kmer frequency. `sum1` means the kmer frequencies are divided by the background kmer frequencies and normalized again to sum to one within each sample/position.
* The kmer size is written into the file name.
* The size of the window is written into the file name.

## Line graph analysis names

The following are parts of file names in the line graph analysis.

* `line_graph`: Line graph of the frequency of ribos along the chromosome.
* `matrix`: Heatmap of the pairwise correlation (or other quantity) of samples along the chromosome. The type of quantity is written into the file name. See the variable `line_graph.SCORE_TYPES` in the `analysis_line_graph.R` script for all possible quantities.
* `outlier_indicators`: Bar graph showing the number of samples that have an outlier at each position on the chromosome. An outlier is defined in one of two ways (see the variable `line_graph.OUTLIER_INDICATOR_FUNCS` in the `analysis_line_graph.R` script for all possible outlier types).

    * `p_value`: A window is an outlier if it is in above a specified empirical quantile. The p-value to determine the quantile is also in the file name (e.g., `0.01` indicates the top 1% are considered outliers).
    * `sd_thresh`: A sample is an outlier if it is more than some standard deviations away from the mean. The number of standard deviations used is also written in the file name.

* `outlier_windows`: Figure showing the regions that are considered outliers in each sample. An outlier is defined in one of two ways (see the variable `line_graph.OUTLIER_INDICATOR_FUNCS` in the `analysis_line_graph.R` script for all possible outlier types).

  * `p_value`: A window is an outlier if it is in above a specified empirical quantile. The p-value to determine the quantile is also in the file name (e.g., `0.01` indicates that the top 1% are considered outliers).
  * `sd_thresh`: A sample is an outlier if it is more than some standard deviations away from the mean. The number of standard deviations used is also written in the file name.

## Neighborhood counts analysis names

The following are parts of file names in the neighborhood counts analysis.

* `multiple`: The figure computes the neighborhood count by weighting each ribo by its pointwise multiplicity (e.g., a ribo with multiplicity 2 is counted twice).
* `single`: The figure computes the neighborhood count by counting each ribo only once, regardless of its multiplicity.
* The size of the window used for computing the neighborhood count is written into the file name.

## Overlap cluster analysis names

The following are parts of file names in the overlap cluster analysis.

* `cluster_n`: Number of clusters at each window size.
* `cluster_n_2`: Number of clusters at each window size, but only for windows with at least 2 ribos.
* `cluster_width_max_2`: Maximum width of clusters at each window size, but only for windows with at least 2 ribos.
* `cluster_width_mean_2`: Mean width of clusters at each window size, but only for windows with at least 2 ribos.
* `cluster_width_median_2`: Median width of clusters at each window size, but only for windows with at least 2 ribos.

## Summary analysis names

The following are parts of file names in the summary analysis. The output from this analysis will be located in a `data_output`-prefixed directory (instead of a `plot`-prefixed directory).

* `reads_per_window`: Contains a table that shows the number of reads-per-window for each sample for different window widths.
* `sample`: Contains a table that shows the cell type, cleavage enzyme, number of ribos, and other information for each sample.
