#!/bin/bash

# NeilRobertson
# Module for indexing bams.  Main requirements are samtools.  Input args = basedirectory


module load igmm/apps/samtools/0.1.19

export LC_ALL=C
export base_directory=$1
export threads=$2

if [ -a $1/.indexedbams ];
then
	echo "Bams already indexed"
else
	echo "Indexing bams"
	# index all remaining bams (4 per sample - all aligned, rm dup removed, multi removed and both removed)
	parallel --progress -j $threads samtools index {1} ::: $1/bam/*.bam

	# move log files to folder
	mv $1/bam/*log $1/logs/

	# remove the sam folder
	#rm -rf $1/sam/
	touch $1/.indexedbams

	echo "Bams indexed"
fi
