#!/bin/bash

# scripts location
heatmap=$1
scripts=$2
output=$3
bg=$4
order=$5

for folder in raw normalized tsv
do
    if ! [ -d $output/heatmap/$folder ]
    then
        mkdir $output/heatmap/$folder
    fi
done

# format
eval $scripts/format.py $output/telomere.tsv -o $output/heatmap/raw/telomere &
eval $scripts/format.py $output/genome_telomere_unit.tsv -o $output/heatmap/raw/genomic &
wait

for ty in both cstrand gstrand
do
    # normalize
    eval $scripts/calc_ef.py $output/heatmap/raw/telomere_mono_${ty}.raw $bg/telomere_mono.raw --name $ty -o $output/heatmap/normalized/telomere_mono_${ty}.norm &
    eval $scripts/calc_ef.py $output/heatmap/raw/telomere_dinuc_d1_nr_${ty}.raw $bg/telomere_dinuc_d1.raw --name $ty -o $output/heatmap/normalized/telomere_dinuc_d1_nr_${ty}.norm &
    eval $scripts/calc_ef.py $output/heatmap/raw/genomic_mono_${ty}.raw $bg/telomere_mono.raw --name $ty -o $output/heatmap/normalized/genomic_mono_${ty}.norm &
    eval $scripts/calc_ef.py $output/heatmap/raw/genomic_dinuc_d1_nr_${ty}.raw $bg/telomere_dinuc_d1.raw --name $ty -o $output/heatmap/normalized/genomic_dinuc_d1_nr_${ty}.norm &
done
wait

# resort
for aa in $(ls $output/heatmap/normalized)
do
    eval $heatmap/resort.py $output/heatmap/normalized/${aa} $order -c 2 -o $output/heatmap/tsv/${aa}.tsv &
done
wait
rename 's/norm.tsv/tsv/' $output/heatmap/tsv/* -f

for ty in cstrand gstrand
do
    # draw
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/telomere_mono_${ty}.tsv -o $output/plots/telomere_mono_${ty}.png --palette RdBu_r --strand ${ty::-6} &
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/telomere_dinuc_d1_nr_${ty}.tsv -o $output/plots/telomere_dinuc_d1_nr_${ty}.png --palette RdBu_r --strand ${ty::-6} &
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/genomic_mono_${ty}.tsv -o $output/plots/genomic_mono_${ty}.png --palette RdBu_r --strand ${ty::-6} &
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/genomic_dinuc_d1_nr_${ty}.tsv -o $output/plots/genomic_dinuc_d1_nr_${ty}.png --palette RdBu_r --strand ${ty::-6} &
done

for ty in both
do
    # draw
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/telomere_mono_${ty}.tsv -o $output/plots/telomere_mono_${ty}.png --palette RdBu_r &
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/telomere_dinuc_d1_nr_${ty}.tsv -o $output/plots/telomere_dinuc_d1_nr_${ty}.png --palette RdBu_r &
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/genomic_mono_${ty}.tsv -o $output/plots/genomic_mono_${ty}.png --palette RdBu_r &
    eval $scripts/draw_heatmap.py $output/heatmap/tsv/genomic_dinuc_d1_nr_${ty}.tsv -o $output/plots/genomic_dinuc_d1_nr_${ty}.png --palette RdBu_r &
done
wait

