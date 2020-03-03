#!/bin/bash

# NeilRobertson

export LC_ALL=C

# move all the fastq files into the sequences directory
if [ -d $1/sequences ];
then
	echo "Sequences data folder already exists - Not moving to sequences folder";
else
	# rename any .txt files into .fastq
	rename 's/.txt/.fastq/g' $1/*.txt*
	
	rename 's/.fq/.fastq/g' $1/*.fq*

	rename 's/_R1_001/_1/g' $1/*.fastq*
	rename 's/_R2_001/_2/g' $1/*.fastq*

	rename 's/_L00//g' $1/*.fastq*

	rename 's/_NoIndex//g' $1/*.fastq*

	rename 's/_R1/_1/g' $1/*.fastq*
	rename 's/_R2/_2/g' $1/*.fastq*

	mkdir $1/sequences
	mv $1/*.fastq $1/sequences/
	mv $1/*.fastq.gz $1/sequences/

	echo "Created restructured/re-labelled sequences folder."
fi

