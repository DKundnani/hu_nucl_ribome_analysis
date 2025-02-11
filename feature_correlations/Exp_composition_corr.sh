#!/bin/bash

#creating reference
source ~/.bash_profile
conda activate r_ecd /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/anno/standardanno/RNAseq
nv

bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -l 499 -r 500 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_500aroundTSS.bed
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -l 999 -r 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_1000aroundTSS.bed
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -b 4000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools flank -s -b 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_4-5aroundTSS.bed
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -b 10000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools flank -s -b 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_30-31aroundTSS.bed

bedtools nuc -fi $hgref -bed sorted_metadata_500aroundTSS.bed -s | tail -n+2  > sorted_metadata_500aroundTSS_withnuc.bed
bedtools nuc -fi $hgref -bed sorted_metadata_1000aroundTSS.bed -s | tail -n+2  > sorted_metadata_1000aroundTSS_withnuc.bed
bedtools nuc -fi $hgref -bed sorted_metadata_4-5aroundTSS.bed -s | tail -n+2  > sorted_metadata_4-5aroundTSS_withnuc.bed
bedtools nuc -fi $hgref -bed sorted_metadata_30-31aroundTSS.bed -s | tail -n+2  > sorted_metadata_30-31aroundTSS_withnuc.bed

#Downstream TSS only
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -l -1 -r 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_1000downTSS.bed
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -l -1 -r 4000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools flank -s -b 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_4-5downTSS.bed
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | grep protein_coding | bedtools slop -s -l -1 -r 9000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools flank -s -b 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_9-10downTSS.bed
bedtools nuc -fi $hgref -bed sorted_metadata_1000downTSS.bed -s | tail -n+2  > sorted_metadata_1000downTSS_withnuc.bed
bedtools nuc -fi $hgref -bed sorted_metadata_4-5downTSS.bed -s | tail -n+2  > sorted_metadata_4-5downTSS_withnuc.bed
bedtools nuc -fi $hgref -bed sorted_metadata_9-10downTSS.bed -s | tail -n+2  > sorted_metadata_9-10downTSS_withnuc.bed


#16_pct_at        17_pct_gc       18_num_A        19_num_C        20_num_G        21_num_T        22_num_N        23_num_oth      24_seq_len

#Separating bed files based on nucleotide
outloc=./poly/
for file in $(ls *.bed); do
  cat $file | awk BEGIN'{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,substr($4,length($4),1)}' > $outloc/poly_$(basename $file)
	for N in A C G T; do
	awk -v n="$N" -F"\t" '{FS=OFS="\t"}{if($7==n) {print $1,$2,$3,$4,$5,$6} }' $outloc/poly_$(basename $file) > $outloc/${N}/$(basename $file)
	done
done


cd /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/
#file='anno/standardanno/RNAseq/sorted_metadata_500aroundTSS_withnuc.bed'; bin=500
file='anno/standardanno/RNAseq/sorted_metadata_1000downTSS_withnuc.bed'; bin=1000
file='anno/standardanno/RNAseq/sorted_metadata_4-5downTSS_withnuc.bed'; bin=4000
#file='anno/standardanno/RNAseq/sorted_metadata_9-10downTSS_withnuc.bed'; bin=30000

subtypecol=1
libmeta='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/heatmaps/libmeta_comp'
genome='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-nucleus-noXY.fa.fai'
ym=12



#Each nucleotide separate
for nuc in A C G T; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
mkdir $outfolder
#bash ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/annotate.sh -r $file -s -c -o $outfolder -b $bedfiles &
sed -i 's/_nucl//g' $outfolder/all_counts.tsv
done


#All rNMPs
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/*.bed')
outfolder=$(echo 'locationHM/raroundTSS'${bin})
mkdir $outfolder
bash ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/annotate.sh -r $file -s -c -o $outfolder -b $bedfiles &

Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/get_celltype_info.R -a $outfolder/annotated_counts.tsv -t locationHM/raroundTSS${bin}/annotated_counts.tsv -c $outfolder/all_counts.tsv -b 16 -f $libmeta -o ${outfolder}/raroundTSS
sed -i 's/_nucl//g' $outfolder/all_counts.tsv
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts.tsv -c $outfolder/all_counts.tsv -g $genome -f $libmeta -t 1 -o ${outfolder}/raroundTSS
<<COMMENT
#Trends from EF
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_EF_avg.tsv -n raroundTSS_regions_WTexp_EF_HEK -c '#D62728' -p -e 13 -r 25 -y $ym -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_EF_avg.tsv -n raroundTSS_regions_KO1exp_EF_HEK -c '#D62728' -p -e 14 -r 25 -y $ym -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_EF_avg.tsv -n raroundTSS_regions_KO2exp_EF_HEK -c '#D62728' -p -e 15 -r 25 -y $ym -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_EF_avg.tsv -n raroundTSS_regions_WTexp_EF_KO1 -c '#E377C2' -p -e 13 -r 26 -y $ym -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_EF_avg.tsv -n raroundTSS_regions_KO1exp_EF_KO1 -c '#E377C2' -p -e 14 -r 26 -y $ym -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_EF_avg.tsv -n raroundTSS_regions_KO2exp_EF_KO1 -c '#E377C2' -p -e 15 -r 26 -y $ym -o locationHM/ribo_exp${bin}/ &
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f raroundTSS_regions_WTexp_EF_ -y 6 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f raroundTSS_regions_KO1exp_EF_ -y 6 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f raroundTSS_regions_KO2exp_EF_ -y 6 -o locationHM/ribo_exp${bin}/
COMMENT

Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_counts_sum.tsv -n raroundTSS_regions_WTexp_raw_HEK -c '#D62728' -p -e 13 -r 25 -y 50 -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_counts_sum.tsv -n raroundTSS_regions_WTexp_raw_KO1 -c '#E377C2' -p -e 13 -r 26 -y 300 -o locationHM/ribo_exp${bin}/ & 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_counts_sum.tsv -n raroundTSS_regions_KO1exp_raw_HEK -c '#D62728' -p -e 14 -r 25 -y 50 -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_counts_sum.tsv -n raroundTSS_regions_KO1exp_raw_KO1 -c '#E377C2' -p -e 14 -r 26 -y 300 -o locationHM/ribo_exp${bin}/ & 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_counts_sum.tsv -n raroundTSS_regions_KO2exp_raw_HEK -c '#D62728' -p -e 15 -r 25 -y 50 -o locationHM/ribo_exp${bin}/ &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/raroundTSS_regions_counts_sum.tsv -n raroundTSS_regions_KO2exp_raw_KO1 -c '#E377C2' -p -e 15 -r 26 -y 300 -o locationHM/ribo_exp${bin}/ &
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f raroundTSS_regions_WTexp_raw_ -y 100 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f raroundTSS_regions_KO1exp_raw_ -y 100 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f raroundTSS_regions_KO2exp_raw_ -y 100 -o locationHM/ribo_exp${bin}/

#Getting normalized percentages
for nuc in A T; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/get_celltype_info.R -a $outfolder/annotated_counts.tsv -t locationHM/raroundTSS${bin}/annotated_counts.tsv -c $outfolder/all_counts.tsv -b 16 -f $libmeta -o ${outfolder}/r${nuc}aroundTSS &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts.tsv -c $outfolder/all_counts.tsv -g $genome -f $libmeta -t 1 -o ${outfolder}/r${nuc}aroundTSS &
done

for nuc in C G; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/get_celltype_info.R -a $outfolder/annotated_counts.tsv -t locationHM/raroundTSS${bin}/annotated_counts.tsv -c $outfolder/all_counts.tsv -b 17 -f $libmeta -o ${outfolder}/r${nuc}aroundTSS &
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts.tsv -c $outfolder/all_counts.tsv -g $genome -f $libmeta -t 1 -o ${outfolder}/r${nuc}aroundTSS &
done


#CD4T = col11, h9 = col12
#Getting trends wrt to HEKExp
for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_raw_HEK -c '#D62728' -p -e 13 -r 25 -y 80 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_raw_KO1 -c '#E377C2' -p -e 13 -r 26 -y 300 -o locationHM/ribo_exp${bin}/  
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_raw_KO2 -c '#9467BD' -p -e 13 -r 27 -y 12 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_perc_HEK -c '#D62728' -p -e 13 -r 25 -y 1 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_perc_KO1 -c '#E377C2' -p -e 13 -r 26 -y 1 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_perc_KO2 -c '#9467BD' -p -e 13 -r 27 -y 12 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_WTexp_norm_HEK -c '#D62728' -p -e 13 -r 25 -y 12 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_WTexp_norm_KO1 -c '#E377C2' -p -e 13 -r 26 -y 12 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_WTexp_norm_KO2 -c '#9467BD' -p -e 13 -r 27 -y 12 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_EF_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_EF_HEK -c '#D62728' -p -e 13 -r 25 -y 12 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_EF_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_EF_KO1 -c '#E377C2' -p -e 13 -r 26 -y 12 -o locationHM/ribo_exp${bin}/ 
done


Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_raw_ -y 250 -o locationHM/ribo_exp${bin}/ &

for bin in 1000 4000 9000; do
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_raw_HEK -y 5 -o locationHM/ribo_exp${bin}/trends 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_raw_KO1 -y 60 -o locationHM/ribo_exp${bin}/trends 

#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_EF_HEK -y 8 -o locationHM/ribo_exp${bin}/trends 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_EF_KO1 -y 8 -o locationHM/ribo_exp${bin}/trends 

#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_perc_HEK -y 1 -o locationHM/ribo_exp${bin}/trends 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/trends_merge.R -i locationHM/ribo_exp${bin}/ -f aroundTSS_regions_WTexp_perc_KO1 -y 1 -o locationHM/ribo_exp${bin}/trends 
done

<<COMMENT

#Getting trends wrt to HEKKO1 Exp
for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_KO1exp_raw_HEK -c '#D62728' -p -e 14 -r 25 -y 100 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_KO1exp_raw_KO1 -c '#E377C2' -p -e 14 -r 26 -y 200 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_KO1exp_raw_KO2 -c '#9467BD' -p -e 14 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_KO1exp_perc_HEK -c '#D62728' -p -e 14 -r 25 -y 1 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_KO1exp_perc_KO1 -c '#E377C2' -p -e 14 -r 26 -y 1 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_KO1exp_perc_KO2 -c '#9467BD' -p -e 14 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_KO1exp_norm_HEK -c '#D62728' -p -e 14 -r 25 -y 8 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_KO1exp_norm_KO1 -c '#E377C2' -p -e 14 -r 26 -y 8 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_KO1exp_norm_KO2 -c '#9467BD' -p -e 14 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_EF_avg.tsv -n r${nuc}aroundTSS_regions_KO1exp_EF_HEK -c '#D62728' -p -e 14 -r 25 -y 12 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_EF_avg.tsv -n r${nuc}aroundTSS_regions_KO1exp_EF_KO1 -c '#E377C2' -p -e 14 -r 26 -y 12 -o locationHM/ribo_exp${bin}/ 
done

for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO1exp_raw_ -y 100 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO1exp_perc_ -y 1 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO1exp_norm_ -y 4 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO1exp_EF_ -y 6 -o locationHM/ribo_exp${bin}/ 
done


#Getting trends wrt to HEKKO2 Exp
for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_KO2exp_raw_HEK -c '#D62728' -p -e 15 -r 25 -y 100 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_KO2exp_raw_KO1 -c '#E377C2' -p -e 15 -r 26 -y 200 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_KO1exp_raw_KO2 -c '#9467BD' -p -e 15 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_KO2exp_perc_HEK -c '#D62728' -p -e 15 -r 25 -y 1 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_KO2exp_perc_KO1 -c '#E377C2' -p -e 15 -r 26 -y 1 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_KO1exp_perc_KO2 -c '#9467BD' -p -e 15 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_KO1exp_norm_HEK -c '#D62728' -p -e 15 -r 25 -y 8 -o locationHM/ribo_exp${bin}/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_KO1exp_norm_KO1 -c '#E377C2' -p -e 15 -r 26 -y 8 -o locationHM/ribo_exp${bin}_2/
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_KO1exp_norm_KO2 -c '#9467BD' -p -e 15 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_EF_avg.tsv -n r${nuc}aroundTSS_regions_KO2exp_EF_HEK -c '#D62728' -p -e 15 -r 25 -y 12 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_EF_avg.tsv -n r${nuc}aroundTSS_regions_KO2exp_EF_KO1 -c '#E377C2' -p -e 15 -r 26 -y 12 -o locationHM/ribo_exp${bin}/ 
done

for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO2exp_raw_ -y 100 -o locationHM/ribo_exp${bin}/ 
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO2exp_perc_ -y 1 -o locationHM/ribo_exp${bin}/ 
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO2exp_norm_ -y 4 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_KO2exp_EF_ -y 6 -o locationHM/ribo_exp${bin}/ 
done

COMMENT
for f in $(ls *_WTexp_*corr.tsv); do
R=$(grep sR $f | cut -f3)
echo -e $f'\t'$R
done

cp ../ribo_exp500/*WTexp*all_line* .
rename 'TSS_' 'TSS500_' *
cp ../ribo_exp1000/*WTexp*all_line* .
rename 'TSS_' 'TSS1000_' *
cp ../ribo_exp4000/*WTexp*all_line* .
rename 'TSS_' 'TSS4000_' *

cp ../ribo_exp500/*WTexp*bins* .
rename 'TSS_' 'TSS500_' *
cp ../ribo_exp1000/*WTexp*bins* .
rename 'TSS_' 'TSS1000_' *
cp ../ribo_exp4000/*WTexp*bins* .
rename 'TSS_' 'TSS4000_' *




