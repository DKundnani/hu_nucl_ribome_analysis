#!/bin/bash

#creating reference
conda activate conda_env
cd /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/anno/standardanno/RNAseq
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | bedtools slop -s -l 499 -r 500 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_500aroundTSS.bed
bedtools flank -s -l 1 -r 0 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i metadata.bed | bedtools slop -s -l 999 -r 1000 -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | bedtools sort -g ~/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38.fa.fai -i stdin | grep -v chrX | grep -v chrY > sorted_metadata_1000aroundTSS.bed
cd /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/
bedtools nuc -fi $hgref -bed anno/standardanno/RNAseq/sorted_metadata_500aroundTSS.bed -s | tail -n+2  > anno/standardanno/RNAseq/sorted_metadata_500aroundTSS_withnuc.bed
bedtools nuc -fi $hgref -bed anno/standardanno/RNAseq/sorted_metadata_1000aroundTSS.bed -s | tail -n+2  > anno/standardanno/RNAseq/sorted_metadata_1000aroundTSS_withnuc.bed

#16_pct_at        17_pct_gc       18_num_A        19_num_C        20_num_G        21_num_T        22_num_N        23_num_oth      24_seq_len

#Separating bed files based on nucleotide
outloc=./poly/
for file in $(ls *.bed); do
  cat $file | awk BEGIN'{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,substr($4,length($4),1)}' > $outloc/poly_$(basename $file)
	for N in A C G T; do
	awk -v n="$N" -F"\t" '{FS=OFS="\t"}{if($7==n) {print $1,$2,$3,$4,$5,$6} }' $outloc/poly_$(basename $file) > $outloc/${N}/$(basename $file)
	done
done

#file='anno/standardanno/RNAseq/sorted_metadata_500aroundTSS_withnuc.bed'; bin=500
file='anno/standardanno/RNAseq/sorted_metadata_1000aroundTSS_withnuc.bed'; bin=1000
subtypecol=1
libmeta='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/heatmaps/libmeta_HEK'
genome='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-nucleus-noXY.fa.fai'

#Each nucleotide separate
for nuc in A C G T; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
mkdir $outfolder
bash ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/annotate.sh -r $file -s -c -o $outfolder -b $bedfiles &
done

#All rNMPs
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/*.bed')
outfolder=$(echo 'locationHM/raroundTSS'${bin})
mkdir $outfolder
bash ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/annotate.sh -r $file -s -c -o $outfolder -b $bedfiles &

#Getting normalized percentages
conda activate r_env
for nuc in A T; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/get_celltype_info.R -a $outfolder/annotated_counts.tsv -t locationHM/raroundTSS${bin}/annotated_counts.tsv -c $outfolder/all_counts.tsv -b 16 -f $libmeta -o ${outfolder}/r${nuc}aroundTSS
done

for nuc in C G; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS'${bin})
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/hu_nucl_ribome_analysis/feature_correlations/get_celltype_info.R -a $outfolder/annotated_counts.tsv -t locationHM/raroundTSS${bin}/annotated_counts.tsv -c $outfolder/all_counts.tsv -b 17 -f $libmeta -o ${outfolder}/r${nuc}aroundTSS
done


#Getting trends wrt to HEKExp
for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS')
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_raw_HEK -c '#D62728' -e 13 -r 25 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_raw_KO1 -c '#E377C2' -e 13 -r 26 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_raw_KO2 -c '#9467BD' -e 13 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_perc_HEK -c '#D62728' -e 13 -r 25 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_perc_KO1 -c '#E377C2' -e 13 -r 26 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_perc_KO2 -c '#9467BD' -e 13 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_WTexp_norm_HEK -c '#D62728' -e 13 -r 25 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_WTexp_norm_KO1 -c '#E377C2' -e 13 -r 26 -y 8 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_norm_perc.tsv -n r${nuc}aroundTSS_regions_WTexp_norm_KO2 -c '#9467BD' -e 13 -r 27 -y 8 -o locationHM/ribo_exp${bin}/
done

for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS')
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_WTexp_raw_ -y 100 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_WTexp_perc_ -y 1 -o locationHM/ribo_exp${bin}/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp${bin}/ -f r${nuc}aroundTSS_regions_WTexp_norm_ -y 4 -o locationHM/ribo_exp${bin}/
done




