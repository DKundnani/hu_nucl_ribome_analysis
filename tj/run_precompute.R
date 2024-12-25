# Initialize
source("channagiri/init.R")
init.init()

# Load the precompute functions
source("channagiri/precompute_hg38_annotations.R")
source("channagiri/precompute_ribos_rds.R")
source("channagiri/precompute_dnaseq_rds.R")
source("channagiri/precompute_samples_test.R")
source("channagiri/precompute_samples_shuffle.R")
source("channagiri/precompute_samples_uniform.R")
source("channagiri/precompute_samples_dnaseq_draw.R")
source("channagiri/precompute_ribos_binned.R")
source("channagiri/precompute_dnaseq_binned.R")
source("channagiri/precompute_weights_binned.R")
source("channagiri/precompute_kmer_binned.R")
source("channagiri/analysis_kmer_heatmaps.R")

# Make the hg38 annotation files
hg38_annotations.do_main(overwrite = FALSE)

# Make the RDS files
ribos_rds.do_main(overwrite = FALSE)
dnaseq_rds.do_main(overwrite = FALSE)

# Make the synthetic samples
samples_test.do_main(overwrite = FALSE)
samples_shuffle.do_main(overwrite = FALSE)
samples_uniform.do_main(overwrite = FALSE)
samples_uniform.do_main_chrAll(overwrite = FALSE)
samples_dnaseq_draw.do_main(overwrite = FALSE)

# Make the sample binned files
ribos_binned.do_main(overwrite = FALSE)
ribos_binned.do_main(test = TRUE, overwrite = FALSE) # For the test data sets
dnaseq_binned.do_main(overwrite = FALSE)

# Make the reference genome binned files
weights_binned.do_main(overwrite = FALSE)
kmer_binned.do_main(overwrite = FALSE)

# Make the heatmap data
kmer_heatmaps.do_precompute(overwrite = FALSE)
kmer_heatmaps.do_precompute(test = TRUE, overwrite = FALSE) # For the test data sets
