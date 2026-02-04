#!/bin/bash

beds=$1
scripts=$2
genome=$3
order=$4
output=$5

for folder in bed_extend fa_extend bed_with_seq
do
    if ! [ -d $output/genome_comparison/$folder ]
    then
        mkdir $output/genome_comparison/$folder
    fi
done

function process {
    name=${1::-4}
    eval $scripts/extend_bed.py $beds/$1 -l 5 | grep -v chrM | grep -v chrX | grep -v chrY | grep -v Done > $output/genome_comparison/bed_extend/$1
    eval bedtools getfasta -fi $genome -bed $output/genome_comparison/bed_extend/$1 -fo $output/genome_comparison/fa_extend/${name}.fa -s -tab
    eval paste $beds/$1 $output/genome_comparison/fa_extend/${name}.fa | cut -f 1,2,3,4,6,8 | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,$4,$6,$5}' > $output/genome_comparison/bed_with_seq/$1
}

# extend
for aa in $(cut -f 2 $order)
do
    file=$aa.bed
    process $file &
done
wait

eval $scripts/rNMP_in_telomere_units.py $output/genome_comparison/bed_with_seq/* $order -o $output/genome_telomere_unit.tsv
echo Done!
