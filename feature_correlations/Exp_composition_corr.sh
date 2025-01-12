#!/bin/bash
outloc=./poly/
for file in $(ls *.bed); do
  cat $file | awk BEGIN'{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,substr($4,length($4),1)}' > $outloc/poly_$(basename $file)
	for N in A C G T; do
	awk -v n="$N" -F"\t" '{FS=OFS="\t"}{if($7==n) {print $1,$2,$3,$4,$5,$6} }' $outloc/poly_$(basename $file) > $outloc/${N}/$(basename $file)
	done
done

file='anno/standardanno/RNAseq/sorted_metadata_500aroundTSS.bed'
subtypecol=1
libmeta='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/heatmaps/libmeta_HEK'
genome='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-nucleus-noXY.fa.fai'

#Each nucleotide separate
for nuc in A C G T; do
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/poly/'${nuc}'/*bed')
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS')
mkdir $outfolder
bash ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/annotate.sh -r $file -s -c -o $outfolder -b $bedfiles &
done

#All rNMPs
bedfiles=$(echo '/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/*.bed')
outfolder=$(echo 'locationHM/raroundTSS')
mkdir $outfolder
bash ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/annotate.sh -r $file -s -c -o $outfolder -b $bedfiles &

#Getting trends

for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS')
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_rawHEK -c '#D62728' -e 13 -r 16 -y 8 -o locationHM/ribo_exp/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_rawKO1 -c '#E377C2' -e 13 -r 17 -y 8 -o locationHM/ribo_exp/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_counts_sum.tsv -n r${nuc}aroundTSS_regions_WTexp_rawKO2 -c '#9467BD' -e 13 -r 18 -y 8 -o locationHM/ribo_exp/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_percHEK -c '#D62728' -e 13 -r 16 -y 8 -o locationHM/ribo_exp/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_percKO1 -c '#E377C2' -e 13 -r 17 -y 8 -o locationHM/ribo_exp/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/pair_bins.R -m ${outfolder}/r${nuc}aroundTSS_regions_perc_avg.tsv -n r${nuc}aroundTSS_regions_WTexp_percKO2 -c '#9467BD' -e 13 -r 18 -y 8 -o locationHM/ribo_exp/
done

for nuc in A C G T; do 
echo $nuc
outfolder=$(echo 'locationHM/r'${nuc}'aroundTSS')
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp/ -f r${nuc}aroundTSS_regions_WTexp_raw -y 20 -o locationHM/ribo_exp/
Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/trends_merge.R -i locationHM/ribo_exp/ -f r${nuc}aroundTSS_regions_WTexp_perc -y 1 -o locationHM/ribo_exp/
done
