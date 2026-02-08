#!usr/bin/env bash

conda activate misc

#conda install bioconda::macs2 #install dependencies if not already installed

if [[ -z $(conda list | grep ucsc-bigwigtobedgraph) ]]; then
  echo "Installing bigWigToBedGraph"
  conda install bioconda::ucsc-bigwigtobedgraph
fi

if [[ -z $(conda list | grep ucsc-bedgraphtobigwig) ]]; then
  echo "Installing bedgraphtobigwig"
  conda install bioconda::ucsc-bedgraphtobigwig
fi

if [[ -z $(conda list | grep macs2) ]]; then
  echo "Installing macs2"
  conda install bioconda::macs2
fi


for file in $(ls *bw); do
echo $file
#bigWigToBedGraph $file $(basename $file .bw).bG # Convert BigWig to bedGraph (requires UCSC tools)
#macs2 callpeak --nomodel --extsize 147 --shift 0 --keep-dup all -g hs -t $(basename $file .bw).bG -q 1e-2 --fe-cutoff 5 -n peaks/$(basename $file .bw) &
macs2 bdgpeakcall -i $(basename $file .bw).bG --o-prefix peaks/$(basename $file .bw)_bgcall --cutoff 10 --min-length 500 &
done
wait

bedtools intersect -a GSM8042853_HEK293T_WT_rep1_pe_redup_bgcall_c5.0_l200_g30_peaks.narrowPeak -b GSM8042854_HEK293T_WT_rep2_pe_redup_bgcall_c5.0_l200_g30_peaks.narrowPeak > bgcall_overlaps.bed
bedtools intersect -a GSM8042853_HEK293T_WT_rep1_pe_redup_peaks.narrowPeak -b GSM8042854_HEK293T_WT_rep2_pe_redup_peaks.narrowPeak > callpeaks_overlap.bed
bedtools intersect -a GSM8042853_HEK293T_WT_rep1_pe_redup_bgcall_c10.0_l500_g30_peaks.narrowPeak -b GSM8042854_HEK293T_WT_rep2_pe_redup_bgcall_c10.0_l500_g30_peaks.narrowPeak > bgcall_overlaps.bed
