#!/bin/bash

# usage
# ./generate_data.sh <fq_folder> <order_list> <output>
fqs=$1
order=$2
output=$3

for folder in raw_reads trimmed trimming_reports logs
do
    if ! [ -e $output/$folder ]
    then
        mkdir $output/$folder
    fi
done

function process {
    # fetch reads
    eval grep 'TTAGGGTTAGGGTTAGGGTTAGGG' $fqs/$1.fq -A 2 -B 1 --no-group-separator > $output/raw_reads/${1}_g.fq &
    eval grep 'CCCTAACCCTAACCCTAACCCTAA' $fqs/$1.fq -A 2 -B 1 --no-group-separator > $output/raw_reads/${1}_c.fq &
    wait
    # # trim reads
    # nohup trim_galore $output/raw_reads/${1}_g.fq -o $output/trimmed > $output/logs/${1}_g.log &
    # nohup trim_galore $output/raw_reads/${1}_c.fq -o $output/trimmed > $output/logs/${1}_c.log &
    # wait 
    # eval "mv $output/trimmed/${1}_*.txt $output/trimming_reports -f"
}

for aa in $(cut -f 2 $order)
do
    process $aa &
done

wait
# rename 's/_trimmed//' $output/trimmed/* -f
