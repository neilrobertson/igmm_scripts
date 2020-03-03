#!/bin/bash

# NeilRobertson
# Module for generating and sorting bams.  Main requirements are samtools.  Input args = basedirectory


module load igmm/apps/samtools/0.1.19

export LC_ALL=C
export base_directory=$1
export threads=$2


# make bams and index them
if [ -d $base_directory/bam ];
then
	echo "BAM data folder already exists - Not creating or indexing BAMs";
else
	echo "Converting to BAM"
	converttobam() {
		# Args:
		# $1 - path to sam file to be converted to bam
		basename=${1%\.*} 				# everything before the last .
		samtools view -bS $1 > $basename.bam 		# make the bam
		samtools sort $basename.bam $basename.sorted ; rm -f $basename.bam ; mv $basename.sorted.bam $basename.bam # sort, remove old and rename
	}
	export -f converttobam
	parallel --progress -j $threads converttobam ::: $1/sam/*.sam

	mkdir $base_directory/bam
	mv $base_directory/sam/*.bam $base_directory/bam/			# move to bam folder

	echo "Completed generating sorted bams from sam."
fi

