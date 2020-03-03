#!/bin/bash

# NeilRobertson
# Module for generating non-normalised counts on eddie.  Main requirements are cufflinks script.  Input args = basedirectory, threads and gtf file

export LC_ALL=C
export base_directory=$1
export threads=$2
export gtf=$3

echo "Using gtf at: ${gtf}"
echo "Using threads: ${threads}"

# Check if alignment completed
if [ -d $1/aligned ];
then
	if [ -d $1/DESeq2 ];
	then
		mkdir $1/DESeq2

		aligned_dirs=$"`ls -d ./aligned/*`"
		for i in $aligned_dirs;
		do
			export rep_dir=${i#"./aligned/"}
			echo $rep_dir

			echo "Beginning to generate non-normalised counts with HTSeq... $rep_dir"

			samtools sort -n -T "$base_directory/aligned/$rep_dir/Aligned.sortedByCoord.out" -@ $threads -o "$base_directory/$rep_dir/DESeq2/$rep_dir.sam" "$base_directory/aligned/$rep_dir/Aligned.sortedByCoord.out.bam"
			samtools view -h -@ $threads "$base_directory/$rep_dir/DESeq2/$rep_dir.sam" > tmp; mv tmp "$base_directory/$rep_dir/DESeq2/$rep_dir.sam"
			htseq-count --mode=intersection-nonempty "$base_directory/$rep_dir/DESeq2/$rep_dir.sam" $gtf > "$base_directory/$rep_dir/DESeq2/$rep_dir.HTSeq.txt"
		done 
	else
		echo "DESeq2 directory already exists. Not generating non-normalised HTSeq Counts"
	fi
else
	echo "Alignment folder not present"
fi
echo "Removing intermediate sam files..."
rm $base_directory/DESeq2/*.sam
echo "Completed HTSeq-Counts generation module"