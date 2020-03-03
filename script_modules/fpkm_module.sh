#!/bin/bash

# NeilRobertson
# Module for generating FPKM values on eddie.  Main requirements are cufflinks script.  Input args = basedirectory, threads and gtf file

export LC_ALL=C
export base_directory=$1
export threads=$2
export gtf=$3

echo "Using gtf at: ${gtf}"
echo "Using threads: ${threads}"

# Check if alignment completed
if [ -d $base_directory/aligned ];
then
	if [ -d $base_directory/FPKM ];
	then
		echo "FPKM directory already exists. Not generating explicit FPKM expression"
	else
		echo "Generating FPKM data..."
		mkdir $base_directory/FPKM

		aligned_dirs=$"`ls -d ./aligned/*`"
		for i in $aligned_dirs;
		do
			export rep_dir=${i#"./aligned/"}
			echo "Beginning to generate FPKM values with Cufflinks... $rep_dir"
			cufflinks_cmd="cufflinks -p $threads -u --max-bundle-frags 1000000 -G $gtf -o $base_directory/FPKM/$rep_dir $base_directory/aligned/$rep_dir/Aligned.sortedByCoord.out.bam"
			echo $cufflinks_cmd
			$cufflinks_cmd
		done 
	fi
else
	echo "Alignment folder not present"
fi

echo "Completed FPKM generation module"