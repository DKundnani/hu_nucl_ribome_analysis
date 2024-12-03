#!usr/bin/env bash

#Wrapper script for Ribosepreference analysis from Penghao
###########################################################
#Get mono, di and tri nucleotide frequencies from BED files overlapping a range file(bed format)
#Usage: source ./Heatmapwrapper.sh
#split_by_subtypes ranges.bed columnusedtosplit

##############################################################################################################


split_by_subtypes () {
        #$1=bed12 file with ranges
        #$2=col number of subtypes
        ranges=$1
        subtypes=$(cut -f${2} ${ranges} | sort | uniq)
        for subtype in ${subtypes}; do grep ${subtype} ${ranges} > ${subtype}_${ranges}; done
}

#intersect multiple bed files with multiple ranges
intersect_multiple () {
        #$1=range files w/ location
        #$2=bedfiles w/ location
        #$3=output folder
        mkdir ${3}
        for range in $(ls ${1}/*bed); do 
        mkdir ${3}/$(basename ${range} .bed)
                for file in $(ls ${2}/*.bed); do
                bedtools intersect -nonamecheck -b ${range} -a ${file} > ${3}/$(basename $range .bed)/$(basename $file)
                done
        done
}

intersect_multiple_ss () {
        mkdir ${3}
        for range in $(ls ${1}/*bed); do 
        mkdir ${3}/$(basename ${range} .bed)
                for file in $(ls ${2}/*.bed); do
                bedtools intersect -s -nonamecheck -b ${range} -a ${file} > ${3}/$(basename $range .bed)/$(basename $file)
                done
        done
}


bg_freq () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #Usage1: bg_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed
        #Usage2: bg_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa ~/p-fstorici3-0/rich_project_bio-storici/AGS/ranges/chrM.bed
        
        
        mkdir -p bg_freq
        bedtools getfasta -fi ${2} -bed ${3} > ./bg_freq/$(basename ${3} .bed).fa
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed).fa -o ./bg_freq/$(basename ${3} .bed).di
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed).di -s chr0 -v -o ./bg_freq/$(basename ${3} .bed).di.freq
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed).fa --mono -o ./bg_freq/$(basename ${3} .bed).mono
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed).mono -s chr0 -v -o ./bg_freq/$(basename ${3} .bed).mono.freq
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed).fa --trinuc -o ./bg_freq/$(basename ${3} .bed).tri
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed).tri -s chr0 -v -o ./bg_freq/$(basename ${3} .bed).tri.freq
}

bg_freq_ss () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #Usage1: bg_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed
        #Usage2: bg_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa ~/p-fstorici3-0/rich_project_bio-storici/AGS/ranges/chrM.bed

        mkdir -p bg_freq
        bedtools getfasta -s -fi ${2} -bed ${3} > ./bg_freq/$(basename ${3} .bed)_same.fa
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed)_same.fa -s -o ./bg_freq/$(basename ${3} .bed)_same.di
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed)_same.di -s chr0 -v -o ./bg_freq/$(basename ${3} .bed)_same.di.freq
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed)_same.fa -s --mono -o ./bg_freq/$(basename ${3} .bed)_same.mono
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed)_same.mono -s chr0 -v -o ./bg_freq/$(basename ${3} .bed)_same.mono.freq
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed)_same.fa -s --trinuc -o ./bg_freq/$(basename ${3} .bed)_same.tri
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed)_same.tri -s chr0 -v -o ./bg_freq/$(basename ${3} .bed)_same.tri.freq
}

bg_freq_os () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #Usage1: bg_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed
        #Usage2: bg_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa ~/p-fstorici3-0/rich_project_bio-storici/AGS/ranges/chrM.bed

        mkdir -p bg_freq
        bedtools getfasta -s -fi ${2} -bed ${3} > ./bg_freq/$(basename ${3} .bed)_opp.fa
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed)_opp.fa -s -o ./bg_freq/$(basename ${3} .bed)_opp.di
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed)_opp.di -s chr0 -v -o ./bg_freq/$(basename ${3} .bed)_opp.di.freq
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed)_opp.fa -s --mono -o ./bg_freq/$(basename ${3} .bed)_opp.mono
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed)_opp.mono -s chr0 -v -o ./bg_freq/$(basename ${3} .bed)_opp.mono.freq
        python3 ${1}/count_background.py ./bg_freq/$(basename ${3} .bed)_opp.fa -s --trinuc -o ./bg_freq/$(basename ${3} .bed)_opp.tri
        python3 ${1}/get_chrom.py ./bg_freq/$(basename ${3} .bed)_opp.tri -s chr0 -v -o ./bg_freq/$(basename ${3} .bed)_opp.tri.freq
}

sample_freq () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/OG_bed/
        mkdir -p sample_freq
        mkdir -p sample_freq/$(basename ${3} .bed)
        for file in $(ls ${4}/*.bed); do
        bedtools intersect -nonamecheck -b ${3} -a ${file} > ./sample_freq/$(basename ${3} .bed)/$(basename ${file})
        #bedtools intersect -s -nonamecheck -b ${3} -a ${file} > ./sample_freq/$(basename ${3} .bed)/$(basename ${file})
        done
        python3 ${1}/count_rNMP.py ${2} ./sample_freq/$(basename ${3} .bed)/*.bed -d -m -t -o ./sample_freq/$(basename ${3} .bed)/
        #python3 ${1}/count_rNMP.py ${2} ./sample_freq/$(basename ${3} .bed)/*.bed --dist 0 -d -m -t -o ./sample_freq/$(basename ${3} .bed)/ #changing distance
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)/*.mono -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_sample.mono
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)/*.dinuc_d1_nr -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_nr
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)/*.dinuc_d1_rn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_rn
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)/*.trinuc_nnr -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_sample.trinuc_nnr
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)/*.trinuc_nrn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_sample.trinuc_nrn
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)/*.trinuc_rnn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_sample.trinuc_rnn
        
}

sample_freq_ss () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/OG_bed/
        mkdir -p sample_freq
        mkdir -p sample_freq/$(basename ${3} .bed)_same
        for file in $(ls ${4}/*.bed); do
        bedtools intersect -s -nonamecheck -b ${3} -a ${file} > ./sample_freq/$(basename ${3} .bed)_same/$(basename ${file})
        done
        python3 ${1}/count_rNMP.py ${2} ./sample_freq/$(basename ${3} .bed)_same/*.bed -d -m -t -o ./sample_freq/$(basename ${3} .bed)_same/
        #python3 ${1}/count_rNMP.py ${2} ./sample_freq/$(basename ${3} .bed)/*.bed --dist 0 -d -m -t -o ./sample_freq/$(basename ${3} .bed)/ #changing distance
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_same/*.mono -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_same.mono
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_same/*.dinuc_d1_nr -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_same.dinuc_d1_nr
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_same/*.dinuc_d1_rn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_same.dinuc_d1_rn
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_same/*.trinuc_nnr -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_same.trinuc_nnr
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_same/*.trinuc_nrn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_same.trinuc_nrn
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_same/*.trinuc_rnn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_same.trinuc_rnn        
}

sample_freq_os () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/OG_bed/
        mkdir -p sample_freq
        mkdir -p sample_freq/$(basename ${3} .bed)_opp
        for file in $(ls ${4}/*.bed); do
        bedtools intersect -S -nonamecheck -b ${3} -a ${file} > ./sample_freq/$(basename ${3} .bed)_opp/$(basename ${file})
        done
        python3 ${1}/count_rNMP.py ${2} ./sample_freq/$(basename ${3} .bed)_opp/*.bed -d -m -t -o ./sample_freq/$(basename ${3} .bed)_opp/
        #python3 ${1}/count_rNMP.py ${2} ./sample_freq/$(basename ${3} .bed)/*.bed --dist 0 -d -m -t -o ./sample_freq/$(basename ${3} .bed)/ #changing distance
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_opp/*.mono -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_opp.mono
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_opp/*.dinuc_d1_nr -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_opp.dinuc_d1_nr
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_opp/*.dinuc_d1_rn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_opp.dinuc_d1_rn
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_opp/*.trinuc_nnr -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_opp.trinuc_nnr
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_opp/*.trinuc_nrn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_opp.trinuc_nrn
        python3 ${1}/get_chrom.py ./sample_freq/$(basename ${3} .bed)_opp/*.trinuc_rnn -s chr0 -v -o ./sample_freq/$(basename ${3} .bed)_opp.trinuc_rnn        
}


norm_freq () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/
        mkdir -p norm_freq
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.mono ./bg_freq/$(basename ${3} .bed).mono.freq -o ./norm_freq/$(basename ${3} .bed)_mono_0 --name $(basename ${3} .bed)
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed).di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_dinuc_nr_4 --name $(basename ${3} .bed)
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed).di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_dinuc_rn_4 --name $(basename ${3} .bed)

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed).di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_dinuc_nr_16 --name $(basename ${3} .bed)
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed).di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_dinuc_rn_16 --name $(basename ${3} .bed)

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.trinuc_nnr ./bg_freq/$(basename ${3} .bed).tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_trinuc_nnr_16 --name $(basename ${3} .bed)
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.trinuc_nrn ./bg_freq/$(basename ${3} .bed).tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_trinuc_nrn_16 --name $(basename ${3} .bed)
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.trinuc_rnn ./bg_freq/$(basename ${3} .bed).tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_trinuc_rnn_16 --name $(basename ${3} .bed)
        
}

norm_freq_ss () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/
        mkdir -p norm_freq
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.mono ./bg_freq/$(basename ${3} .bed)_same.mono.freq -o ./norm_freq/$(basename ${3} .bed)_same_mono_0 --name $(basename ${3} .bed)_same
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed)_same.di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_same_dinuc_nr_4 --name $(basename ${3} .bed)_same
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed)_same.di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_same_dinuc_rn_4 --name $(basename ${3} .bed)_same

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed)_same.di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_same_dinuc_nr_16 --name $(basename ${3} .bed)_same
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed)_same.di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_same_dinuc_rn_16 --name $(basename ${3} .bed)_same

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.trinuc_nnr ./bg_freq/$(basename ${3} .bed)_same.tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_same_trinuc_nnr_16 --name $(basename ${3} .bed)_same
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.trinuc_nrn ./bg_freq/$(basename ${3} .bed)_same.tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_same_trinuc_nrn_16 --name $(basename ${3} .bed)_same
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_same.trinuc_rnn ./bg_freq/$(basename ${3} .bed)_same.tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_same_trinuc_rnn_16 --name $(basename ${3} .bed)_same
        
}

norm_freq_os () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/
        mkdir -p norm_freq
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.mono ./bg_freq/$(basename ${3} .bed)_opp.mono.freq -o ./norm_freq/$(basename ${3} .bed)_opp_mono_0 --name $(basename ${3} .bed)_opp
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed)_opp.di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_opp_dinuc_nr_4 --name $(basename ${3} .bed)_opp
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed)_opp.di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_opp_dinuc_rn_4 --name $(basename ${3} .bed)_opp

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed)_opp.di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_opp_dinuc_nr_16 --name $(basename ${3} .bed)_opp
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed)_opp.di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_opp_dinuc_rn_16 --name $(basename ${3} .bed)_opp

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.trinuc_nnr ./bg_freq/$(basename ${3} .bed)_opp.tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_opp_trinuc_nnr_16 --name $(basename ${3} .bed)_opp
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.trinuc_nrn ./bg_freq/$(basename ${3} .bed)_opp.tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_opp_trinuc_nrn_16 --name $(basename ${3} .bed)_opp
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_opp.trinuc_rnn ./bg_freq/$(basename ${3} .bed)_opp.tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_opp_trinuc_rnn_16 --name $(basename ${3} .bed)_opp
        
}

norm_freq_noname () {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #Usage1: sample_freq ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/
        mkdir -p norm_freq
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.mono ./bg_freq/$(basename ${3} .bed).mono.freq -o ./norm_freq/$(basename ${3} .bed)_mono_0 

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed).di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_dinuc_nr_4 
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed).di.freq --group_len 4 -o ./norm_freq/$(basename ${3} .bed)_dinuc_rn_4 

        #python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_nr ./bg_freq/$(basename ${3} .bed).di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_dinuc_nr_16 
        #python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.dinuc_d1_rn ./bg_freq/$(basename ${3} .bed).di.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_dinuc_rn_16 

        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.trinuc_nnr ./bg_freq/$(basename ${3} .bed).tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_trinuc_nnr_16 
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.trinuc_nrn ./bg_freq/$(basename ${3} .bed).tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_trinuc_nrn_16 
        python3 ${1}/normalize.py ./sample_freq/$(basename ${3} .bed)_sample.trinuc_rnn ./bg_freq/$(basename ${3} .bed).tri.freq --group_len 16 -o ./norm_freq/$(basename ${3} .bed)_trinuc_rnn_16 
        
}

resort_plot() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        #Usage1: resort_plot ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/ 
        for file in $(ls ./norm_freq/$(basename ${3} .bed)*[0-9]); do 
        python3 ${1}/resort.py ${file} ${5} -o ./norm_freq/sorted_$(basename $file)
        done

        mkdir -p plots 
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*mono*); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed).mono.freq --background_chrom $(basename ${3} .bed) -o ./plots/$(basename $file) --palette RdBu_r --group_size 4
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*dinuc*4); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed).di.freq --background_chrom $(basename ${3} .bed) -o ./plots/$(basename $file) --palette RdBu_r --group_size 4
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*dinuc*16); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed).di.freq --background_chrom $(basename ${3} .bed) -o ./plots/$(basename $file) --palette RdBu_r --group_size 16
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*trinuc*16); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed).tri.freq --background_chrom $(basename ${3} .bed) -o ./plots/$(basename $file) --palette RdBu_r --group_size 16
        done

}

resort_plot_ss() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        #Usage1: resort_plot ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/ 
        for file in $(ls ./norm_freq/$(basename ${3} .bed)*same*[0-9]); do 
        python3 ${1}/resort.py ${file} ${5} -o ./norm_freq/sorted_$(basename $file)
        done

        mkdir -p plots 
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*same*mono*); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_same.mono.freq --background_chrom $(basename ${3} .bed)_same -o ./plots/$(basename $file)_same --palette RdBu_r --group_size 4
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*same*dinuc*4); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_same.di.freq --background_chrom $(basename ${3} .bed)_same -o ./plots/$(basename $file)_same --palette RdBu_r --group_size 4
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*same*dinuc*16); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_same.di.freq --background_chrom $(basename ${3} .bed)_same -o ./plots/$(basename $file)_same --palette RdBu_r --group_size 16
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*same*trinuc*16); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_same.tri.freq --background_chrom $(basename ${3} .bed)_same -o ./plots/$(basename $file)_same --palette RdBu_r --group_size 16
        done

}

resort_plot_os() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        #Usage1: resort_plot ~/p-fstorici3-0/rich_project_bio-storici/bin/RibosePreferenceAnalysis/ ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/anno/subtypes/hg38_cpg_islands.bed ~/p-fstorici3-0/rich_project_bio-storici/HEKnH9/bed/ 
        for file in $(ls ./norm_freq/$(basename ${3} .bed)*opp*[0-9]); do 
        python3 ${1}/resort.py ${file} ${5} -o ./norm_freq/sorted_$(basename $file)
        done

        mkdir -p plots 
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*opp*mono*); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_opp.mono.freq --background_chrom $(basename ${3} .bed)_opp -o ./plots/$(basename $file)_opp --palette RdBu_r --group_size 4
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*opp*dinuc*4); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_opp.di.freq --background_chrom $(basename ${3} .bed)_opp -o ./plots/$(basename $file)_opp --palette RdBu_r --group_size 4
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*opp*dinuc*16); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_opp.di.freq --background_chrom $(basename ${3} .bed)_opp -o ./plots/$(basename $file)_opp --palette RdBu_r --group_size 16
        done
        for file in $(ls ./norm_freq/sorted_$(basename ${3} .bed)*opp*trinuc*16); do
        python3 ${1}/draw_heatmap.py ${file} -b ./bg_freq/$(basename ${3} .bed)_opp.tri.freq --background_chrom $(basename ${3} .bed)_opp -o ./plots/$(basename $file)_opp --palette RdBu_r --group_size 16
        done

}

mww() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        mkdir mww
        SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_mono_0 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_mono_0_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_dinuc_nr_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_dinuc_nr_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_dinuc_rn_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_dinuc_rn_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_nnr_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_trinuc_nnr_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_nrn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_trinuc_nrn_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_rnn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_trinuc_rnn_16_mww_pval.tsv
}

mww_same() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        mkdir mww
        SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_same_mono_0 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_same_mono_0_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_same_dinuc_nr_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_same_dinuc_nr_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_same_dinuc_rn_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_same_dinuc_rn_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_same_trinuc_nnr_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_same_trinuc_nnr_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_same_trinuc_nrn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_same_trinuc_nrn_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_same_trinuc_rnn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_same_trinuc_rnn_16_mww_pval.tsv
}

mww_opp() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        #$6=liborder
        #$7=column of labels in liborder
        #$8=column of groups in liborder
        mkdir mww
        SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_opp_mono_0 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_opp_mono_0_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_opp_dinuc_nr_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_opp_dinuc_nr_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_opp_dinuc_rn_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww/$(basename ${3} .bed)_opp_dinuc_rn_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_opp_trinuc_nnr_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_opp_trinuc_nnr_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_opp_trinuc_nrn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_opp_trinuc_nrn_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww.R -f ./norm_freq/sorted_$(basename ${3} .bed)_opp_trinuc_rnn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww/$(basename ${3} .bed)_opp_trinuc_rnn_16_mww_pval.tsv
}


ttest() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of ranges
        #$4=location of bed files
        #$5=order file
        mkdir ttest
        SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        Rscript $SCRIPT_DIR/ttest.R -f ./norm_freq/sorted_$(basename ${3} .bed)_mono_0 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./ttest/$(basename ${3} .bed)_mono_0_mww_pval.tsv
        Rscript $SCRIPT_DIR/ttest.R -f ./norm_freq/sorted_$(basename ${3} .bed)_dinuc_nr_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./ttest/$(basename ${3} .bed)_dinuc_nr_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/ttest.R -f ./norm_freq/sorted_$(basename ${3} .bed)_dinuc_rn_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./ttest/$(basename ${3} .bed)_dinuc_rn_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/ttest.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_nnr_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./ttest/$(basename ${3} .bed)_trinuc_nnr_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/ttest.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_nrn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./ttest/$(basename ${3} .bed)_trinuc_nrn_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/ttest.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_rnn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./ttest/$(basename ${3} .bed)_trinuc_rnn_16_mww_pval.tsv
}

mww_diff() {
        #$1=Location of Scripts/git repo
        #$2=Ref fasta
        #$3=bed12 file of range 1
        #$4=bed12 file of range 2
        #$5=order file
        #$6=liborder
        #$7=column of labels in liborder
        #$8=column of groups in liborder
        mkdir mww_diff
        SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        Rscript $SCRIPT_DIR/mww_diff.R -f ./norm_freq/sorted_$(basename ${3} .bed)_mono_0 -F ./norm_freq/sorted_$(basename ${4} .bed)_mono_0 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww_diff/$(basename ${3} .bed)_mono_0_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww_diff.R -f ./norm_freq/sorted_$(basename ${3} .bed)_dinuc_nr_4 -F ./norm_freq/sorted_$(basename ${4} .bed)_dinuc_nr_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww_diff/$(basename ${3} .bed)_dinuc_nr_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww_diff.R -f ./norm_freq/sorted_$(basename ${3} .bed)_dinuc_rn_4 -F ./norm_freq/sorted_$(basename ${4} .bed)_dinuc_rn_4 -m ${6} -n ${7} -b 0.25 -g ${8} -o ./mww_diff/$(basename ${3} .bed)_dinuc_rn_4_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww_diff.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_nnr_16 -F ./norm_freq/sorted_$(basename ${4} .bed)_trinuc_nnr_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww_diff/$(basename ${3} .bed)_trinuc_nnr_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww_diff.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_nnr_16 -F ./norm_freq/sorted_$(basename ${4} .bed)_trinuc_nnr_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww_diff/$(basename ${3} .bed)_trinuc_nnr_16_mww_pval.tsv
        Rscript $SCRIPT_DIR/mww_diff.R -f ./norm_freq/sorted_$(basename ${3} .bed)_trinuc_rnn_16 -F ./norm_freq/sorted_$(basename ${4} .bed)_trinuc_rnn_16 -m ${6} -n ${7} -b 0.0625 -g ${8} -o ./mww_diff/$(basename ${3} .bed)_trinuc_rnn_16_mww_pval.tsv
}