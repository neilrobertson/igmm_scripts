#!/bin/bash

# NeilRobertson
# Module for a hisat2 alignment on eddie.   Input args = basedirectory

export LC_ALL=C
export base_directory=$1
export threads=$2
export index=$3
export quals="auto"

echo "Using index at: ${index}"
echo "Using threads: ${threads}"

export PYTHONPATH="/mnt/tchandra-lab/Neil/repository/workspace/Baselib/:$PYTHONPATH"
export fastqdetectpath="/mnt/tchandra-lab/Neil/repository/workspace/Tools/fastq/fastq-detect.py"

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

		if [ $quals == "auto" ];
        	then
            		qual="`python $fastqdetectpath --fastq $i`"
            		# solexa quals is actually phred64 here
            		if [ $qual == "solexa" ];
            		then
                		qual="phred64"
            		fi
        	else
            		qual=$quals
        	fi

		hisatlog="$i.hisat2.log"

		date > $hisatlog
		echo "***"
		echo "Perfoming alignment on file $alignmentNumber. Filename: $i"
		echo "HISAT2 log filename: $hisatlog"

		hisat2 --version >> $hisatlog

		one=$i
		two=`echo $i | sed "s/_1/_2/"`

		basename=${one##*/}
		name=${basename%_*}

		outfolder=$1"/aligned/"$name"/"
		mkdir $outfolder
		echo "Outfolder: $outfolder"

		if [ -f $two ]
		then
			hisatcmd="hisat2 --dta-cufflinks -x ${index} --${qual} -p ${threads} -q -S $outfolder/accepted_hits.sam -1 $one -2 $two"
		else
			hisatcmd="hisat2 --dta-cufflinks -x ${index} --${qual} -p ${threads} -q -S $outfolder/accepted_hits.sam -U $one"
		fi

		echo $hisatcmd >> $hisatlog
		echo "***"
		echo $hisatcmd
		echo "***"
		$hisatcmd &>> $hisatlog
	done

	# move tophat logs to logs folder
	if [ -d $1/logs ];
	then
		echo "Logs folder exists"
	else
		mkdir $1/logs
	fi
	mv $1/trimmed-seq/*log $1/logs/

	echo "Completed HISAT2 alignment module."
fi
