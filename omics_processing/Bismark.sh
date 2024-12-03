#Initializing variables
conda activate conda-env #activate conda environment with bismark tool 
ref='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/hg38/filtered_hg38' #reference fasta file
aligned='/storage/home/hcoda1/5/dkundnani3/scratch/Bismark/bams' #folder for alignment files
trimmed='/storage/home/hcoda1/5/dkundnani3/scratch/Bismark/trimmed' #Location of trimmed reads
lib='SampleN' #It should be the sample name of the library
in='path/to/${lib}/' #Location of the input fastq files

###################################################
# Start pipeline
###################################################
###### Prepare genome
bismark_genome_preparation --bowtie2 $ref
###### Trimming
R1=$in/${lib}_R1.fq.gz
R2=$in/${lib}_R2.fq.gz
FQ1=${trimmed}/$(basename $R1 .fq)_val_1.fq
FQ2=${trimmed}/$(basename $R2 .fq)_val_2.fq
####Taking small sample ~1mil reads of the fastq files (optional)

#cd $in
#seqtk sample -s100 ${SLURM_ARRAY_TASK_ID}_R1.fastq.gz 10000000 > ${SLURM_ARRAY_TASK_ID}_sub1.fq 
#seqtk sample -s100 ${SLURM_ARRAY_TASK_ID}_R2.fastq.gz 10000000 > ${SLURM_ARRAY_TASK_ID}_sub2.fq 
#R1=$in/${lib}_sub1.fq
#R2=$in/${lib}_sub2.fq
#FQ1=${trimmed}/$(basename $R1 .fq)_val_1.fq
#FQ2=${trimmed}/$(basename $R2 .fq)_val_2.fq

trim_galore -q 20 --cores 8 --fastqc -o $trimmed --clip_r1 8 --clip_r2 20 --three_prime_clip_r1 8 --three_prime_clip_r2 8 --paired ${R1} ${R2}

#Unzip trimmed fq files
#cd $trimmed
gunzip ${FQ1}.gz 
gunzip ${FQ2}.gz

###################################################
# Run Bismark 
###################################################
cd ${aligned}
bismark --bowtie2 --bam --multicore 6 $ref -1 $FQ1 -2 $FQ2
#bismark --non_directional --bam --multicore 6 --bowtie2 -p 2 $ref -1 $FQ1 -2 $FQ2 #optional check
deduplicate_bismark --bam ${lib}_R1_val_1_bismark_bt2_pe.bam
bismark_methylation_extractor -p --bedgraph --zero_based --counts --buffer_size 20G --cytosine_report --genome_folder $ref --multicore 10 --gzip --ignore 5 --ignore_r2 5 ${lib}_R1_val_1_bismark_bt2_pe.deduplicated.bam
#bismark2report #optionanl check
