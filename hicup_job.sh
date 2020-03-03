#!/bin/sh

# Set replicate run-time parameters:
export replicate_name=`basename "$PWD"`  #Grabs the current folder and assumes sequences are labelled similarly

export species="GRCh38"
export digestion="/mnt/tchandra-lab/Neil/genomes/GRCh38/restriction_enzymes/Digest_GRCh38_HindIII_None_11-47-50_24-10-2018.txt"
export threads="8"

####################################################
################## Begin actual work ###############
####################################################

construct_sequences_folder.sh .

fastqc_module.sh . $threads

trimseqs_module.sh . $threads auto

# Run HiCUP alignment
if [ -d ./hicup/ ];
then
	echo "HiCUP directory exists - not aligning HiC with HiCUP.";
else
	# Designate the hicup tool path:

	mkdir ./hicup/
	
	hicup --zip \
	--bowtie2 /usr/bin/bowtie2 \
	--index /mnt/tchandra-lab/Neil/genomes/${species}/INDEXES/${species} \
	--digest $digestion \
	--format Sanger \
	--longest 800 \
	--shortest 150 \
	--threads $threads \
	--outdir ./hicup/ \
	./trimmed-seq/${replicate_name}_1.fastq.gz ./trimmed-seq/${replicate_name}_2.fastq.gz

	echo "HiCUP alignment completed."
fi

# Run HOMER HiCzummary conversion.
if [ -s ./hicup/${replicate_name}_1_2.hicup.bam.homer.gz ];
then
	echo "HiCUP to HOMER format conversion exists. Skipping."
else
	echo "Beginning homer output generation..."
	/mnt/tchandra-lab/Neil/repository/3rdparty/hicup_v0.6.1/Conversion/hicup2homer --zip ./hicup/${replicate_name}_1_2.hicup.bam
	echo "Homer output completed."

fi


# Run HOMER conversion.
if [ -d ./homer/ ];
then
	echo "homer folder exists. Skipping.";
else
	echo "Beginning generation of homer tag directory (makeTagDirectory)..."
	makeTagDirectory homer -format HiCsummary ./hicup/${replicate_name}_1_2.hicup.bam.homer.gz -tbp 1
	echo "Homer tag directory completed."

	echo "Generating .hic file......"
	# Add juicer jar to path
	tagDir2hicFile.pl ./homer -juicer auto -genome hg38 -p ${threads}
	echo ".hic file generated!"
fi

