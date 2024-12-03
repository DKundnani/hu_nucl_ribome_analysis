#!usr/bin/env bash

###################MODIFY AS PER THE BATCH RUNS ###################################
scripts='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/' #Locaiton for the Preference analysis scripts/ local github repository
ref='human-genome/non-altchromosomes.fa' #Fasta file for the reference genome
range='human-genome/nucl.bed' #bed file of ranges in the nucleus of choice of ranges
bed='ribodataset/bedfolder/' #folder of bed files with extension .bed
order='ribodataset/order' #Sequence of used bed file names and what you would like to append to it. Sequence should reflect the order in which you would like to have the files plotted!

###################DO NOT MODIFY BELOW THIS LINE ###################################
bg_freq $scripts $ref $range 
sample_freq $scripts $ref $range $bed
norm_freq $scripts $ref $range $bed
resort_plot $scripts $ref $range $bed $order