#!/bin/bash

# usage
# ./generate_data.sh <trimmed_fq_folder> <SAM folder> <order_list> <output>
fqs=$1
sam=$2
order=$3
output=$4

for folder in telomere telomeric_genome
do
    if ! [ -e $output/$folder ]
    then
        mkdir $output/$folder
    fi
done

function process {
    # fetch telomere reads
    eval cat $fqs/${1}_R1.fq $fqs/${1}_R2.fq \
        | grep 'TTAGGGTTAGGGTTAGGGTTAGGG' \
        > $output/telomere/${1}_g.txt &
    eval cat $fqs/${1}_R1.fq $fqs/${1}_R2.fq \
        | grep 'CCCTAACCCTAACCCTAACCCTAA' \
        > $output/telomere/${1}_c.txt &
    # fetch telomeric region reads in genome
    eval samtools view -F 2308 $sam/${1}.sam \
        | cut -f 10 \
        | grep 'TTAGGG' \
        > $output/telomeric_genome/${1}_g.txt &
    eval samtools view -F 2308 $sam/${1}.sam \
        | cut -f 10 \
        | grep 'CCCTAA' \
        > $output/telomeric_genome/${1}_c.txt &
    wait
}

# Process individual reads
for aa in $(cut -f 2 $order)
do
    process $aa &
done
wait

# summarize
current=$(pwd)
cd $output/telomere
wc -l *.txt | grep -v total | sed 's/^ *//;s/ /\t/;s/.txt//;s/_/\t/' \
    | awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' > ../DNA_telomere_count.tsv
cd ../telomeric_genome
wc -l *.txt | grep -v total | sed 's/^ *//;s/ /\t/;s/.txt//;s/_/\t/' \
    | awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' > ../DNA_telomeric_genome_count.tsv
cd $current
