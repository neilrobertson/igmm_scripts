#!/bin/bash

# NeilRobertson
# Module for a STAR alignment on eddie.  Main requirements are trim-galore/cutadapt and fastq-detect script.  Input args = basedirectory

echo $3
export LC_ALL=C
export base_directory=$1
export threads=$2
export species=$3
export gtf=$4

export genome_path="/mnt/tchandra-lab/Neil/genomes/"${species}"/CELLRANGER_INDEX/"
echo "Using index at: ${genome_path}"
echo "Using threads: ${threads}"

# do the alignment
if [ -d $1/aligned ];
then
	echo "Alignment already done - Not aligning";
else
	mkdir $1/aligned
	export alignmentCount=0
	for i in $1/trimmed-seq/*_1*.fastq*
	do
		alignmentCount=$((alignmentCount+1))
		alignmentNumber=$( printf '%03d' $alignmentCount )
		outfile=$i$mapsuffix

		starlog="$i.star.log"

		date > $starlog
		echo "***"
		echo "Perfoming alignment on file $alignmentNumber. Filename: $i"
		echo "STAR log filename: $tophatlog"

		STAR --version >> $starlog

		one=$i
		two=`echo $i | sed "s/_1/_2/"`

		# GTF file example (hg19) : ftp://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz
		# remember to rewrite it to add "chr" to chromosome names if not present already
		
		# cd /mnt/bam01/publicdata/hg19
		# wget ftp://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz
		# gunzip Homo_sapiens.GRCh37.73.gtf.gz
		# cat Homo_sapiens.GRCh37.73.gtf | sed 's/\(.*\)/chr\1/' > Homo_sapiens.GRCh37.73.fixedchr.gtf

		# outfolder = ./aligned/<name>
		basename=${one##*/}
		name=${basename%_*}

		outfolder=$1"/aligned/"$name"/"
		mkdir $outfolder
		echo "Outfolder: $outfolder"

		if [ -f $two ]
		then
			starcmd="STAR --runThreadN "${threads}" \
			--genomeDir "${genome_path}" \
			--outSAMtype BAM SortedByCoordinate \
			--readFilesCommand zcat \
			--quantMode GeneCounts \
			--outSAMstrandField intronMotif \
			--outFilterMultimapNmax 1 \
			--readFilesIn $one $two \
			--outFileNamePrefix $outfolder"
		else
			starcmd="STAR --runThreadN ${threads} \
			--genomeDir ${genome_path} \
			--outSAMtype BAM SortedByCoordinate \
			--readFilesCommand zcat \
			--quantMode GeneCounts \
			--outSAMstrandField intronMotif \
			--outFilterMultimapNmax 1 \
			--readFilesIn $one \
			--outFileNamePrefix $outfolder"
		fi

		echo $starcmd >> $starlog
		echo "***"
		echo $starcmd
		echo "***"
		$starcmd &>> $starlog
	done

	# move tophat logs to logs folder
	if [ -d $1/logs ];
	then
		echo "Logs folder exists"
	else
		mkdir $1/logs
	fi
	mv $1/trimmed-seq/*log $1/logs/

	echo "Completed STAR alignment module."
fi
