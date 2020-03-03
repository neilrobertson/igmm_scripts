
#!/bin/bash

# NeilRobertson
# Module for trimming sequences on eddie.  Main requirements are trim-galore/cutadapt and fastq-detect script.  Input args = basedirectory

export LC_ALL=C
export base_directory=$1
export threads=$2
export quals=$3

echo "Use threads: $threads"


export PYTHONPATH="/mnt/tchandra-lab/Neil/repository/workspace/Baselib/:$PYTHONPATH"
export fastqdetectpath="/mnt/tchandra-lab/Neil/repository/workspace/Tools/fastq/fastq-detect.py"

if [ -d $base_directory/trimmed-seq ];
then
	echo "Trimmed sequence exists - Not trimming";
else
	mkdir $base_directory/trimmed-seq

	echo "Trimming sequences"

	trimseq() {
		# Args:
		# $1 - path to file to be trimmed
		# $2 - path for output

		echo "Trimming $quals"

		if [ $quals == "auto" ];
		then
			qual="`python $fastqdetectpath --fastq $1`"

			# solexa quals is actually phred64 here
			if [ $qual == "solexa" ];
			then
				qual="phred64"
			fi
		else
			qual=$quals
		fi
		
		one=$1
		two=`echo $1 | sed "s/_1/_2/"`

		if [ -f $two ]
		then
			/mnt/tchandra-lab/Neil/repository/3rdparty/TrimGalore-0.5.0/trim_galore --quality 30 --trim-n --$qual --paired --suppress_warn --retain_unpaired $one $two --output_dir $2
		else
			/mnt/tchandra-lab/Neil/repository/3rdparty/TrimGalore-0.5.0/trim_galore --quality 30 --trim-n --$qual --suppress_warn $one --output_dir $2
		fi
	}
	export -f trimseq

	parallel --progress -j $threads trimseq ::: $base_directory/sequences/*_1*.fastq* ::: $base_directory/trimmed-seq

	if [ -d $base_directory/trimmed-seq/unpaired ];
	then
		echo "Not making unpaired folder"
	else
		mkdir $base_directory/trimmed-seq/unpaired
	fi

	mv $base_directory/trimmed-seq/*unpaired*fq* $base_directory/trimmed-seq/unpaired/

	rename 's/_1_val_1.fq.gz/_1.fastq.gz/g' $base_directory/trimmed-seq/*_1_val_1.fq.gz
	rename 's/_2_val_2.fq.gz/_2.fastq.gz/g' $base_directory/trimmed-seq/*_2_val_2.fq.gz

	rename _trimmed.fq .fastq $base_directory/trimmed-seq/*fq*

	if [ -d $base_directory/logs ];
	then
		echo "Logs folder exists"
	else
		mkdir $base_directory/logs
	fi
	mv $base_directory/trimmed-seq/*trimming_report.txt $base_directory/logs/
	
	echo "Read trimming module completed."
fi
