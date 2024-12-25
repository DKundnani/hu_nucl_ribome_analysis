# Initialize
source("channagiri/init.R")
init.init()

# Source the analysis functions
source("channagiri/analysis_summary.R")
source("channagiri/analysis_kmer_heatmaps.R")
source("channagiri/analysis_kmer_correlation.R")
source("channagiri/analysis_kmer_hist_offset.R")
source("channagiri/analysis_distribution.R")
source("channagiri/analysis_line_graph.R")
source("channagiri/analysis_neighborhood_counts.R")
source("channagiri/analysis_between_ribo_dist.R")
source("channagiri/analysis_dnaseq.R")
source("channagiri/analysis_overlap_clusters.R")

# Run the analyses
summary.do_main(test = settings.TEST)
kmer_heatmaps.do_main(test = settings.TEST, overwrite = TRUE)
kmer_correlation.do_main(test = settings.TEST)
# kmer_hist_offset.do_main(test = settings.TEST) # commented out because it takes too long to run and not being used
distribution.do_main(test = settings.TEST)
line_graph.do_main(test = settings.TEST)
neighborhood_counts.do_main(test = settings.TEST)
between_ribo_dist.do_main(test = settings.TEST)
dnaseq_analysis.do_main(test = settings.TEST)
overlap_clusters.do_main(test = settings.TEST)
