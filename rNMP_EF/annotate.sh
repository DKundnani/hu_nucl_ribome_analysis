#! usr/bin/bash

get_input () {
	
	while getopts "r:o:fchisdb" option
	do 
		case ${option} in
			r) reference=${OPTARG}; ;;
			b) b=1; ;;
			o) output=${OPTARG}; ;;
			c) counts=1; ;;
			f) coverage=1; ;;
			s) strand=1; ;;
			d) divide=1; ;;
			i) intensity; ;;
			h) echo "usage: bash annotate.sh -r reference_annotation.bed -c -o output -b bedfile(s)
					[-r <file with reference annotations, make sure the file is sorted>] 
					[-b <bed file name, wild card enabled for multiple bed files>] 
					[-c <calculate counts in the referenced annotations>] 
					[-f <calculate fraction of coverage in the reference annotations>]
					[-s <annotations are strand specific>]
					[-r <file with reference annotations>] 
					[-o <output file name, will be stored in the annotations folder created in the current directory>]
					[-d <annotated files will be separated based on nucleus and mitochondria>]
					[-h] help option, gives usage information 
				";;
			*) 

				
		esac
	done
	shift $(( OPTIND - 1 ))
	bed=$@
	echo 
	
	#flag check
	if [[ $reference == "" ]]; then
	echo "Please provide a reference"
	exit 1
	fi

	if [[ $b == "" ]]; then
	echo "Please provide bed files "
	exit 1
	fi

	if [[ $counts == "" ]] && [[ $coverage == "" ]]; then
	echo "Please specify what you want to calculate: counts of coveragae fraction"
	exit 1
	fi
  
}


initiate () {
	#check for reference annotation and bed files 
	if [[ -f "$reference" ]]; then echo "$reference found"; else echo "annotation file not found in the specified location"; fi
	echo 
	#for file in "$files"
	for file in ${bed} 
	do
	if [[ -f $file ]];then echo "$file found"; else "Bed files not found in the current locations"; fi 
	done
	echo 
	#create output directory
	if [[ -d $output ]]; then echo "Removing already existing output directory"; rm -r $output;mkdir $output; [[ -d "annotations" ]] && echo "New output directory created"  ; else mkdir $output; [[ -d $output ]] && echo "output directory created"; fi
	echo
	#check if the reference annotation format is compatible for strandness
	if [[ "$strand" == 1 ]]; then col6=$(cut -f6 "$reference" | head -1); if [[ ${col6} =~ "+" ]] || [[ "${col6}" =~ "-" ]]; then echo "Reference annotation format looks good"; else echo "Reference format not compatible. Please have bed6 file with strand information in the 6th column"; fi; else echo "No strand information needed"; fi
  	if [[ $divide == 1 ]]; then echo "Diving by organelle"; else echo "Not diving by organelle"; fi
}


annotating() {
	##### Getting counts #####
	if [[ $counts == 1 ]]; then
	echo
	echo "Annotating"
	echo
	bedtools annotate -counts -i $reference -files $bed > $output/annotated_counts.tsv
	echo "Annotating counts for both strands complete"
		if [[ $strand == 1 ]]; then
		#same strand
		bedtools annotate -counts -s -i $reference -files $bed > $output/annotated_counts_same.tsv
		echo "Annotating counts on same strand complete"
		#opposite strand
		bedtools annotate -counts -S -i $reference -files $bed > $output/annotated_counts_opp.tsv
		echo "Annotating counts on opposite strand complete"
		fi 
	fi 
	##### Getting fraction coverage #####
	if [[ $coverage == 1 ]]; then
	bedtools annotate -i $reference -files $bed > $output/annotated_frac.tsv
	echo "Annotating coverage/fraction for both strands complete"
		if [[ $strand == 1 ]]; then
		#same strand
		bedtools annotate -s -i $reference -files $bed > $output/annotated_frac_same.tsv
		echo "Annotating coverage/fraction on same strand complete"
		#opposite strand
		bedtools annotate -S -i $reference -files $bed > $output/annotated_frac_opp.tsv
		echo "Annotating coverage/fraction on opposite strand complete"
		fi 
	fi 
	##### Getting intensity: Counts/(Coverage fraction * Length)
	#awk '{$11 = $9 / $10}1' file | sort -nr -k 11 | column -t
}

split() {
	if [[ $divide == 1 ]]; then
	mkdir $output/chrM
  echo "Created chrM directory in the output directory"
  mkdir $output/nucl
  echo "Created nucl directory in the output directory"
  echo 
 
  for file in $(ls $output/*.tsv)
	do
	grep "^chrM" ${file} > $output/chrM/chrM_$(basename $file)
	echo "Annotations for chrM as stored in new "chrM" directory"
	#Compatible with Human and Sacceromyces
	grep -v "^chrM" ${file} | grep -v "^2micron" >  $output/nucl/nucl_$(basename $file) #works for both yeast and human
	echo "Annotations for nucleus as stored in new "nucl" directory"
	done
  echo "Finished splitting annotation files as per organelles"
	fi
	
  
}

allcounts() {
  touch $output/all_counts.tsv
	for file in ${bed} 
	do
  filename=$(basename $file)
  count=$(cat $file | wc -l)
  echo $filename $count >> $output/all_counts.tsv
	done
  
  if [[ $divide == 1 ]]; then
  touch $output/chrM/chrM_all_counts.tsv
  touch $output/nucl/nucl_all_counts.tsv
  for file in ${bed} 
	do
  filename=$(basename $file)
	chr_count=$(grep "^chrM" ${file} | wc -l)
  echo $filename $chr_count>> $output/chrM/chrM_all_counts.tsv
  nucl_count=$(grep -v "^chrM" ${file} | grep -v "^2micron" | wc -l)
  echo $filename $nucl_count>> $output/nucl/nucl_all_counts.tsv
	done
  fi
 
}

main() {
	get_input "$@"
	initiate 
	annotating 
  split
  allcounts
}

# Calling the main function
main "$@"