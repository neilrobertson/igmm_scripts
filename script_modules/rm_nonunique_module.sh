#!/bin/bash

# NeilRobertson
# Module for removing non-unique alignments.  Main requirements are early samtools.  Input args = basedirectory, cores-require

export LC_ALL=C
export base_directory=$1
export threads=$2

module load python/2.7.10
module load igmm/apps/samtools/0.1.19


if [ -a $base_directory/.nonuniqremoved ];
then
	echo "Multimapped reads already removed"
else
	echo "Removing multimapped reads"
	# make bams with multimapped reads removed (will also remove from previous rmdup removed)
	removemultimapped() {
		# Args:
		# $1 - path to bam file to remove multimapped reads from
		samtools view -h $1 | grep -Fv "XS:i:" | samtools view -b -S -o $1.uniq -
	}
	export -f removemultimapped
	parallel --progress -j $threads removemultimapped ::: $base_directory/bam/*.bam
	rename .bam.uniq .uniq.bam $base_directory/bam/*.uniq
	touch $base_directory/.nonuniqremoved

	echo "Removed multi-mapped reads (.uniq.bam)"
fi

