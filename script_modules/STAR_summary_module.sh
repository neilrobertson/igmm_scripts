#!/bin/bash

# NeilRobertson
# Module for an RNA-seq alignment via STAR on eddie.  Main requirements are python script.  Input args = basedirectory


export LC_ALL=C
export base_directory=$1
export threads=$2

export script_path="/home/nerobert/bin/script_modules"

echo "Using base directory at: ${base_directory}"

# Check if alignment completed
if [ -d $base_directory/aligned ];
then

	echo "Generating STAR summary data..."
	mkdir alignment_stats
	python $script_path/STAR_get_alignment_stats.py --input "${base_directory}/aligned" --output "${base_directory}/alignment_stats/"

	mkdir counts
	python $script_path/STAR_combine_counts.py --input "${base_directory}/aligned" --output "${base_directory}/counts/" --strand "1"
else
	echo "Alignment folder not present"
fi

# Check if FPKM generation completed
if [ -d $base_directory/FPKM ];
then

	echo "Summarizing and generating FPKM matrix..."
	python $script_path/FPKM_combine_counts.py --FPKM_directory "${base_directory}/FPKM"
else
	echo "FPKM folder not present"
fi

echo "Completed STAR alignment summary module"