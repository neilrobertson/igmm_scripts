#!/bin/bash

export dir=$1
export threads=$2
export gtf=$3

cd aligned
echo $PWD

conda activate

htseq_make() {
	export i="./${1}"
	samtools sort -n -T ${i}/accepted_hits -o ${i}/accepted_hits.sorted.sam ${i}/accepted_hits.sam
	samtools view -h ${i}/accepted_hits.sorted.sam > ${i}.tmp; mv ${i}.tmp ${i}/accepted_hits.sorted.sam
	htseq-count --mode=intersection-nonempty ${i}/accepted_hits.sorted.sam ${gtf} > ${i}/accepted_hits.htseq-counts.txt
}

export -f htseq_make

parallel --progress -j ${threads} htseq_make ::: *

cd ..
conda deactivate

cd aligned
processbams() {
	export i="./${1}"
	samtools view -b -S -o ${i}/accepted_hits.bam ${i}/accepted_hits.sam
	samtools sort -o ${i}/accepted_hits.sorted.bam ${i}/accepted_hits.bam ; rm -f ${i}/accepted_hits.bam ; mv ${i}/accepted_hits.sorted.bam ${i}/accepted_hits.bam # sort, remove old and rename
	samtools index ${i}/accepted_hits.bam
}

export -f processbams

parallel --progress -j ${threads} processbams ::: *
cd ..
