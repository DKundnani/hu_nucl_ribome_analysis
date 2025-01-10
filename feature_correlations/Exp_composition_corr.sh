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



Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts.tsv -c $outfolder/all_counts.tsv -g $genome -f $libmeta -t $subtypecol -o ${outfolder}/$(basename $file)
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts_same.tsv -c $outfolder/all_counts.tsv -f $libmeta -t $subtypecol  -o ${outfolder}/same_$(basename $file)
#Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts_opp.tsv -c $outfolder/all_counts.tsv -f $libmeta -t $subtypecol  -o ${outfolder}/opp_$(basename $file)
