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

	echo "Generating Kallisto summary data..."
	mkdir counts
	python $script_path/Kallisto_combine_counts.py --input "${base_directory}/aligned" --output "${base_directory}/counts/"
else
	echo "Alignment folder not present"
fi

echo "Completed Kallisto alignment summary module"