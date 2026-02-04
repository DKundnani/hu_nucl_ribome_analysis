#!/bin/bash

# usage
script=$1
fqs=$2
order=$3
output=$4
beds=$5
genome=$6
heatmap=$7

# create folders
for folder in reads plots genome_comparison heatmap
do
    if ! [ -d $output/$folder ]
    then
        mkdir $output/$folder
    fi
done

# generate reads
eval $script/generate_data.sh $fqs $order $output/reads
echo "Telomere Data generated!"

# generate barchart
eval "$script/barchart/loc_freq.py $output/reads/raw_reads/*.fq $order -o $output/telomere.tsv --umi_length 0"
eval $script/barchart/draw_freq.py $output/telomere.tsv -o $output/plots/telomere_barchart.png
echo 'Telomere bar chart generated!'

# generate rNMP in genome
eval $script/genome_comparison/generate_comparison.sh $beds $script/genome_comparison \
    $genome $order $output
eval $script/barchart/draw_freq.py $output/genome_telomere_unit.tsv -o $output/plots/genome_barchart.png
echo 'genome bar chart generated!'

# generate heatmap
eval $script/heatmap/draw.sh $heatmap $script/heatmap $output $script/heatmap/bg $order
echo 'Heatmap generated!'

wait
