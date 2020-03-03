#!/bin/bash

# NeilRobertson
# Module for trimming sequences on eddie.  Main requirements are fastQC and fastq-detect script.  Input args = basedirectory and number of cores

export LC_ALL=C
export base_directory=$1
export threads=$2

export contaminants_list="/mnt/tchandra-lab/Neil/repository/3rdparty/fastqc/contaminants_list.txt"

# if the quality folder already exists we can skip, otherwise call FastQC
if [ -d $base_directory/quality ];
then
	echo "Quality data folder already exists - Not quality checking";
else
	echo "Quality reports. Cores: ${threads}"
	mkdir $base_directory/quality
	/mnt/tchandra-lab/Neil/repository/3rdparty/fastqc_v0.11.8/FastQC/fastqc --nogroup --threads $threads --contaminants $contaminants_list -f fastq $base_directory/sequences/*.fastq $base_directory/sequences/*.fastq.gz
	mv $base_directory/sequences/*fastqc* $base_directory/quality/

	echo "Completed FastQC quality check. Resuts in ./quality"
fi
