#!/bin/bash

# Draw heatmaps and barplots
# Usage ./heatmap_barplot.sh <RibosePreferenceAnalysis location> <nuc genome> <bed folder> <library order> <output folder> 
RibosePreferenceAnalysis=$1
genome=$2
bed_folder=$3
order=$4
output_folder=$5

# make folders
if ! [ -d $output_folder ]
then
    mkdir $output_folder
fi

for folder in individual raw normalized tsv plots bg
do
    if ! [ -d $output_folder/$folder ]
    then
        mkdir $output_folder/$folder
    fi
done

# background
eval $RibosePreferenceAnalysis/count_background.py $genome --mono -o $output_folder/bg/nuc_mono.raw1 &
for dis in $(seq 1 5) 100
do
    eval $RibosePreferenceAnalysis/count_background.py $genome -d ${dis} -o $output_folder/bg/nuc_dinuc_d${dis}.raw1 &
done
eval $RibosePreferenceAnalysis/count_background.py $genome --trinuc -o $output_folder/bg/nuc_trinuc.raw1 &
wait

ls $output_folder/bg/*.raw1 | xargs -I aa $RibosePreferenceAnalysis/get_chrom.py aa -v -s chrM chrX chrY -a --name nuc -o aa.raw
rename 's/raw1.raw/raw/' $output_folder/bg/* -f
rm $output_folder/bg/*.raw1

# count individual
for file in $(ls $bed_folder)
do
    eval $RibosePreferenceAnalysis/count_rNMP.py $genome $bed_folder/$file -m -d --dist 1 2 3 4 5 100 -t -o $output_folder/individual/nuc &
done
wait
rename 's/nuc_//' $output_folder/individual/* -f

# Get all files
eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*mono -v -s chrM chrX chrY -o $output_folder/raw/nuc_mono.raw &
for dis in 1 2 3 4 5 100
do
    eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*d${dis}_nr -v -s chrM chrX chrY -o $output_folder/raw/nuc_dinuc_d${dis}_nr.raw &
    eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*d${dis}_rn -v -s chrM chrX chrY -o  $output_folder/raw/nuc_dinuc_d${dis}_rn.raw &
done
for pat in nnr nrn rnn 
do
    eval $RibosePreferenceAnalysis/get_chrom.py $output_folder/individual/*${pat} -v -s chrM chrX chrY -o $output_folder/raw/nuc_trinuc_${pat}.raw &
done
wait

# normalize
eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/nuc_mono.raw $output_folder/bg/nuc_mono.raw --name nuc -o $output_folder/normalized/nuc_mono.norm &
for dis in 1 2 3 4 5 100
do
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/nuc_dinuc_d${dis}_nr.raw $output_folder/bg/nuc_dinuc_d${dis}.raw --name nuc --group_len 4 -o $output_folder/normalized/nuc_dinuc_d${dis}_nr.norm &
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/nuc_dinuc_d${dis}_rn.raw $output_folder/bg/nuc_dinuc_d${dis}.raw --name nuc --group_len 4 -o $output_folder/normalized/nuc_dinuc_d${dis}_rn.norm &
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/nuc_dinuc_d${dis}_nr.raw $output_folder/bg/nuc_dinuc_d${dis}.raw --name nuc -o $output_folder/normalized/nuc_dinuc_d${dis}_nr_16.norm &
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/nuc_dinuc_d${dis}_rn.raw $output_folder/bg/nuc_dinuc_d${dis}.raw --name nuc -o $output_folder/normalized/nuc_dinuc_d${dis}_rn_16.norm &
done
for pat in nnr nrn rnn
do
    eval $RibosePreferenceAnalysis/normalize.py $output_folder/raw/nuc_trinuc_${pat}.raw $output_folder/bg/nuc_trinuc.raw --name nuc --group_len 16 -o $output_folder/normalized/nuc_trinuc_${pat}.norm &
done
wait

# resort
ls $output_folder/normalized | xargs -I aa -P 8 $RibosePreferenceAnalysis/resort.py $output_folder/normalized/aa $order -c 1 -o $output_folder/tsv/aa.tsv
rename 's/norm.tsv/tsv/' $output_folder/tsv/* -f

# draw
eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/nuc_mono.tsv -b $output_folder/bg/nuc_mono.raw -o $output_folder/plots/nuc_mono.png --palette RdBu_r &
for dis in 1 2 3 4 5 100
do
    for ty in nr rn
    do
        eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/nuc_dinuc_d${dis}_${ty}.tsv -b $output_folder/bg/nuc_dinuc_d${dis}.raw -o $output_folder/plots/nuc_dinuc_d${dis}_${ty}.png --palette RdBu_r &
        eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/nuc_dinuc_d${dis}_${ty}_16.tsv -b $output_folder/bg/nuc_dinuc_d${dis}.raw -o $output_folder/plots/nuc_dinuc_d${dis}_${ty}_16.png --palette RdBu_r &
    done
done
for pat in nnr nrn rnn
do
    eval $RibosePreferenceAnalysis/draw_heatmap.py $output_folder/tsv/nuc_trinuc_${pat}.tsv -b $output_folder/bg/nuc_trinuc.raw -o $output_folder/plots/nuc_trinuc_${pat}.png --no_annot --palette RdBu_r &
done
wait


# draw barplot
eval $RibosePreferenceAnalysis/generate_bar_plot.py $output_folder/tsv/nuc_mono.tsv -o $output_folder/plots/nuc_barplot_normalized.png&
eval $RibosePreferenceAnalysis/sum1.py $output_folder/raw/nuc_mono.raw -o $output_folder/normalized/nuc_sum1.norm
eval $RibosePreferenceAnalysis/resort.py $output_folder/normalized/nuc_sum1.norm $order -c 1 -o $output_folder/tsv/nuc_mono_sum1.tsv
eval $RibosePreferenceAnalysis/generate_bar_plot.py $output_folder/tsv/nuc_mono_sum1.tsv -o $output_folder/plots/nuc_barplot_raw.png &
