#!/bin/bash

# NeilRobertson
# Module for trimming sequences on eddie.  Main requirements are trim-galore/cutadapt and fastq-detect script.  Input args = basedirectory

module load igmm/apps/bowtie/2.3.1
module load python/2.7.10
module load igmm/apps/pigz/2.3.3


export LC_ALL=C
export base_directory=$1
export threads=$2
export species=$3
export quals=$4

export BOWTIE_PATH="/exports/igmm/eddie/tchandra-lab/Neil/genomes/$species/INDEXES/"
export BOWTIE2_PATH="/exports/igmm/eddie/tchandra-lab/Neil/genomes/$species/INDEXES/"
export PYTHONPATH="/exports/igmm/eddie/tchandra-lab/Neil/repository/workspace/Baselib/:$PYTHONPATH"

export fastqdetectpath="/exports/igmm/eddie/tchandra-lab/Neil/repository/workspace/Tools/fastq/fastq-detect.py"

# do the alignment
if [ -d $base_directory/sam ];
then
	echo "SAM data folder already exists - Not aligning";
else
	if [ -d $base_directory/bam ];
	then
		echo "BAM data folder already exists - Not aligning";
	else
		mkdir $base_directory/sam

		for i in $base_directory/trimmed-seq/*_1*.fastq*
		do
			outfile=$i$mapsuffix

			bowtielog="$i.bowtie.log"

			date > $bowtielog
			echo $bowtielog

			bowtie2 --version >> $bowtielog

			if [ $quals == "auto" ];
			then
				# note we use the original fastq file for quality detection as we might change the score distribution by trimming
				# i.e. only retain high scoring reads all of which have scores in the phred64 range even though it's actually phred33
				basename=${i##*/}
				qual="`python $fastqdetectpath --fastq $base_directory/sequences/$basename`"

				# solexa quals needs "-quals" added onto the end
				if [ $qual == "solexa" ];
				then
					qual="solexa-quals"
				fi
			else
				qual=$quals
			fi

			echo "Using quality: $qual"

			one=$i
			two=`echo $i | sed "s/_1/_2/"`

			if [ -f $two ]
			then
				bowtiecmd="bowtie2 --$qual -p $threads -x $species -1 $one -2 $two -S $outfile"
			else
				bowtiecmd="bowtie2 --$qual -p $threads -x $species -U $one -S $outfile"
			fi
			
			echo $bowtiecmd >> $bowtielog
			echo "***"
			echo $bowtiecmd
			echo "***"
			$bowtiecmd &>> $bowtielog

			mv $outfile $base_directory/sam/
		done

		# move bowtie logs to logs folder
		if [ -d $base_directory/logs ];
		then
			echo "Logs folder exists"
		else
			mkdir $base_directory/logs
		fi
		mv $base_directory/trimmed-seq/*log $base_directory/logs/
	fi
	echo "Completed bowtie2 alignment stage.  Outputs in sam directory"
fi