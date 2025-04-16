#finding overlap of CpG islands with other annotations
for file in $(ls ../../subtypes/*.bed); do
overlap=$(grep -v chrX cpg_away_TSS_1kb_around.bed |  grep -v chrY | bedtools intersect -wao -a stdin -b $file| grep -e '\schr' | cut -f4 | sort | uniq | wc -l)
tots=$(grep -v chrX cpg_away_TSS_1kb_around.bed |  grep -v chrY | wc -l)
echo -e $file"\t"$(echo $overlap/$tots*100 | bc -l)
done

#getting coordinates for CpG islands overlapping with other annotations
for file in $(ls ../../subtypes/*.bed); do
grep -v chrX cpg_away_TSS_center.bed |  grep -v chrY | bedtools intersect -wao -a stdin -b $file| grep -e '\schr' | cut -f1-6 | bedtools sort -g $hggenome -i stdin | uniq > $(basename $file)
done

#Finding overlap of overlaps (which CpG overlapping annotations have overlaps?)
for file in $(ls *.bed); do
echo "~~~~~~~~~~~~~~~~"
echo $file
echo "~~~~~~~~~~~~~~~~"
for f in $(ls *.bed); do
overlap=$(bedtools intersect -wao -a $file -b $f | grep -e '\schr' | cut -f4 | sort | uniq | wc -l)
tots=$(cat $file | wc -l)
echo -e $f"\t"$(echo $overlap/$tots*100 | bc -l)
done
done


#studying overlap combining ccre and gene annotations

grep -v chrX cpg_away_TSS_center.bed |  grep -v chrY | bedtools intersect -wao -a stdin -b ../../preprocessed/ccre-hg38.bed | grep -e '\schr' | cut -f17 | sort | uniq -c

grep -v chrX cpg_away_TSS_center.bed |  grep -v chrY | bedtools intersect -wao -a stdin -b ../../preprocessed/ccre-hg38.bed | grep -ve '\schr' | cut -f1-6 | bedtools intersect -wao -a stdin -b ../../proteincoding_expandedexons.bed_short | grep -e '\schr' | cut -f13 | sort | uniq -c

grep -v chrX cpg_away_TSS_center.bed |  grep -v chrY | bedtools intersect -wao -a stdin -b ../../preprocessed/ccre-hg38.bed | grep -ve '\schr' | cut -f1-6 | bedtools intersect -wao -a stdin -b ../../proteincoding_expandedexons.bed_short | grep -ve '\schr' | cut -f1-6 | bedtools intersect -wao -a stdin -b ../../standardanno/combined_genic_anno_noXY.bed | cut -f13,16 | sort | uniq -c


#Location percent heatmaps rNMPs in 1 kb around CpG centers aways from TSS
#cpg away from TSS
cd /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/cpgawayTSS
for file in $(ls ../*.bed); do grep -v chrX ../../../../anno/standardanno/cpg_awayTSS/cpg_away_TSS_center.bed | grep -v chrY | bedtools slop -b 1000 -i stdin -g $hggenome | bedtools intersect -wa -a $file -b stdin > $(basename $file) & done

conda activate r_env #RE
libmeta='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/files_short'
bedfiles='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/Hu_analysis/subnfiltbed/nucl/noXY/cpgawayTSS/*.bed' #cpgawayTSS rNMPs
genome='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38-nucleus-noXY.fa.fai'

locationHM () {
    bash ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/annotate.sh -r ${file} -c -o $outfolder -b $bedfiles
    Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/TAVIR/Subtype_ratio.R -a $outfolder/annotated_counts.tsv -c $outfolder/all_counts.tsv -g $genome -f $libmeta -t $subtypecol -o ${outfolder}/$(basename $file)
    Rscript ~/p-fstorici3-0/rich_project_bio-storici/bin/GIT/TAVIR/Subtype_ration_plots.R -o $libmeta -f $outfolder/$(basename $file)_subtype_percent.tsv -y 20
}

cd /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis
file='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/Hu_analysis/anno/proteincoding_expandedintronexons.bed'
outfolder='locationHM/cpgawaTSS_pc_expandedintronexons'
subtypecol=7
locationHM

file='anno/standardanno/combined_ccre_annno.bed'
outfolder='locationHM/cpgawayTSS_ccre'
subtypecol=11
locationHM
