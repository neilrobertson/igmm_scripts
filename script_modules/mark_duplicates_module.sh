#!/bin/bash

# NeilRobertson
# Module for removing non-unique alignments.  Main requirements are early samtools.  Input args = basedirectory, cores-require

export LC_ALL=C
export base_directory=$1
export threads=$2

module load python/2.7.10
module load igmm/apps/samtools/0.1.19
module load igmm/apps/jdk/1.8.0_66
module load igmm/apps/picard/1.139



if [ -a $base_directory/.dupsremoved ];
then
	echo "Dups already removed"
else
	echo "Removing dups"
	# make a bam with duplicates removed (PCR or optical)
	removedups() {
		# Args:
		# $1 - path to bam file to remove duplicates reads from
		java -Xms256m -Xmx4g -jar /exports/igmm/software/pkg/el7/apps/picard/1.139/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=$1 OUTPUT=$1.rmdup METRICS_FILE=$1.dupmetrics.log ASSUME_SORTED=true &> $1.markdup.log
	}
	export -f removedups
	parallel --progress -j $threads removedups ::: $base_directory/bam/*.bam
	rename .bam.rmdup .rmdup.bam $base_directory/bam/*.rmdup
	touch $base_directory/.dupsremoved

	echo "Scripts has marked duplicates with picard tools"
fi
