#!/usr/bin/env bash

#download the initiation index big wig files from GSE189668 and name them as ii_rep1.bw ii_rep2.bw ...
#modify the names of files in the following script as needed. 

#Getting bedgraph files
bigWigToBedGraph ii_rep1.bw ii_rep1.bg
bigWigToBedGraph ii_rep2.bw ii_rep2.bg
bigWigToBedGraph ii_rep3.bw ii_rep3.bg

#Getting bedgraph initiation indexes on one file
bedtools unionbedg -i ii_rep1.bg ii_rep2.bg ii_rep3.bg -o sum > merged.bg

#Avergae ii for each file
awk 'BEGIN { OFS = "\t" }{ print $1, $2, $3, ($4+$5+$6) / 3 }' merged.bg > averaged.bed

#Getting those regions with common overlap
bedtools intersect -wa -c -a averaged.bed -b ii_rep1.bg ii_rep2.bg ii_rep3.bg | grep -e '\s3' > common_regions.bed

# Cluster the regions in the top 10% BED file for 250,000 bases
bedtools cluster -i common_regions.bed -d 100000 > clustered_common.bed

# Report the highest value in each cluster
awk 'BEGIN { OFS = "\t" } { if ($4 > max[$6]) { max[$6] = $4; line[$6] = $0 } } END { for (i in line) print line[i] }' clustered_common.bed > highest_in_clusters.bed

## Extract the top 10% lines
awk '$4 > 0' highest_in_clusters.bed | sort -k4,4nr > sorted_highestincluster.bed

total_lines=$(wc -l < highest_in_clusters.bed)
top_10_percent_lines=$((total_lines / 4))

head -n $top_10_percent_lines highest_in_clusters.bed | sort -k1,1 -k2,2n > top_10_highestincluster.bed

## Use bedtools to get 1kb 5 kb and 10 kb around top 10% ii's
for bin in 1000 5000 10000
do
    bedtools slop -i top_10_highestincluster.bed -g $genome -l $bin -r 0 -s | cut -f1-4,6 | awk 'BEGIN { OFS = "\t" } { print $0, "+" }' > top_10_ii_${bin}bp_up.bed
    bedtools slop -i top_10_highestincluster.bed -g $genome -l 0 -r $bin -s | cut -f1-4,6 | awk 'BEGIN { OFS = "\t" } { print $0, "+" }' > top_10_ii_${bin}bp_down.bed
done

# Get heatmaps
scripts='../RibosePreferenceAnalysis/' #Locaiton for the Preference analysis scripts/ local github repository
ref=../../genome/pu-seq/filtered_hg38-nucleus-noXY.fa #Fasta file for the reference genome
#ref='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/reference/sacCer3/sacCer3.fa' #Fasta file for the reference genome
bed=../bed/
order=order_pole
source ../RPA-wrapper/Heatmapwrapper.sh
conda activate misc

for range in $(ls ../range/*); do
#bg_freq_ss $scripts $ref $range ; bg_freq_os $scripts $ref $range
sample_freq_ss $scripts $ref $range $bed ; sample_freq_os $scripts $ref $range $bed
#norm_freq_ss $scripts $ref $range $bed ; norm_freq_os $scripts $ref $range $bed
#resort_plot_ss $scripts $ref $range $bed $order ; resort_plot_os $scripts $ref $range $bed $order
#rm sample_freq/$(basename $range .bed)/*bed #Cleanup
#mww_same $scripts $ref $range $bed $order $libmeta 6 4
#mww_opp $scripts $ref $range $bed $order $libmeta 6 4
#mww_diff $scripts $ref $(basename ${range} .bed)_same $(basename ${range} .bed)_opp $order $libmeta 6 4
done


