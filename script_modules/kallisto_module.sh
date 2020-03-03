#!/bin/bash

# NeilRobertson
# Module for a Kallisto alignment on eddie.  Main requirements are trim-galore/cutadapt and fastq-detect script.  Input args = basedirectory

echo $3
export LC_ALL=C
export base_directory=$1
export threads=$2
export species=$3
export index=$4
export gtf=$5

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

		kallistolog="$i.kallisto.log"

		date > $kallistolog
		echo "***"
		echo "Perfoming alignment on file $alignmentNumber. Filename: $i"
		echo "Kallisto log filename: $kallistolog"

		kallisto version >> $kallistolog

		one=$i
		two=`echo $i | sed "s/_1/_2/"`

		basename=${one##*/}
		name=${basename%_*}

		outfolder=$1"/aligned/"$name"/"
		mkdir $outfolder
		echo "Outfolder: $outfolder"

		if [ -f $two ]
		then
			kallistocmd="kallisto quant --index=${index} \
			--output-dir=${outfolder} \
			--bias \
			--gtf ${gtf} \
			--threads=${threads} \
			$one $two"
		else
			kallistocmd="kallisto quant --index=${index} \
			--output-dir=${outfolder} \
			--bias \
			--gtf ${gtf} \
			--threads=${threads} \
			$one"
		fi

		echo $kallistocmd >> $kallistolog
		echo "***"
		echo $kallistocmd
		echo "***"
		$kallistocmd &>> $kallistolog
	done

	# move tophat logs to logs folder
	if [ -d $1/logs ];
	then
		echo "Logs folder exists"
	else
		mkdir $1/logs
	fi
	mv $1/trimmed-seq/*log $1/logs/

	echo "Completed Kallisto alignment module."
fi
