#!/usr/bin/env bash

#download the initiation index big wig files from GSE189668 and name them as ii_rep1.bw ii_rep2.bw ...
#modify the names of files in the following script as needed. 

bigWigToBedGraph ii_rep1.bw ii_rep1.bg
bigWigToBedGraph ii_rep2.bw ii_rep2.bg
bigWigToBedGraph ii_rep3.bw ii_rep3.bg


bedtools unionbedg -i ii_rep1.bg ii_rep2.bg ii_rep3.bg -o sum > merged.bg

awk 'BEGIN { OFS = "\t" }{ print $1, $2, $3, ($4+$5+$6) / 3 }' merged.bg > averaged.bed

bedtools intersect -wa -c -a averaged.bed -b ii_rep1.bg ii_rep2.bg ii_rep3.bg | grep -e '\s3' > common_regions.bed