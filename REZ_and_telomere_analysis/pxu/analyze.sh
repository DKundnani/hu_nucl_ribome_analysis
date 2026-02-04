#!/bin/bash

# usage
# scripts
script='~/bio-storici/scripts/rNMP_human_Nuclear_analysis/'
heatmap='~/bio-storici/scripts/RibosePreferenceAnalysis/'
fastq_folder='/storage/home/hcoda1/0/pxu64/scratch/raw_reads'
bed_folder='~/bio-storici/shared/bed/filtered/human/'
telomere_order="$script/telomere/order.tsv"
libinfo="~/bio-storici/scripts/rNMP_human_Nuclear_analysis/libinfo.tsv"
genome="~/bio-storici/genome/hg38.fa"

output='results_updated'

# create folders
if ! [ -d $output ]
then
    mkdir $output
fi

for folder in telomere plots transcript basic logs
do
    if ! [ -e $output/$folder ]
    then
        mkdir $output/$folder
    fi
done

# # Count rNMPs in each library
# if [ -e $output/nuc_count.tsv ]
# then
#     rm $output/nuc_count.tsv
# fi
# curr=$(pwd)
# eval cd $bed_folder
# for aa in $(ls)
# do
#     grep -v chrM $aa | grep -v chrX | grep -v chrY | wc -l | paste <(echo $aa) - | sed 's/.bed//' >> $curr/$output/nuc_count.tsv
# done
# cd $curr

# Analyze telomere
eval $script/telomere/analyze.sh $script/telomere $fastq_folder $telomere_order $output/telomere \
    $bed_folder $genome $heatmap 2>&1 > $output/logs/telomere.log &

# # Analyze pyramid
# eval $script/transcript/analyze_rNMP_in_transcript.sh $script/transcript $bed_folder $genome $libinfo $output/nuc_count.tsv $output/transcript 2>&1 > $output/logs/transcript.log &
# eval $script/transcript/analyze_rNMP_pyramid_longer_bins.sh \
#     $script/transcript $bed_folder $genome $libinfo $output/nuc_count.tsv \
#     $output/transcript_longer_bins 2>&1 \
#     > $output/logs/transcript_longer_bins.log &

# # basic study
# eval $script/analyze_heatmap_barplot.sh $heatmap $genome $bed_folder $libinfo $output/basic 2>&1 > $output/logs/basic.log &

wait

# # move figures
# mv $output/telomere/plots $output/plots/telomere -f 
# mv $output/transcript/plots $output/plots/transcript -f 
# mv $output/transcript/pyramid_analysis/plots $output/plots/pyramid -f 
# mv $output/transcript_longer_bins/pyramid_analysis/plots/* $output/plots/pyramid -f 
# mv $output/basic/plots $output/plots/basic -f

echo "Done!"