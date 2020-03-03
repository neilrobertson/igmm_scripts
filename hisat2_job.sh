#!/bin/bash

export PATH="/home/nerobert/bin/script_modules:$PATH"

#Set project specific params
export genome_build="GRCh38"
export gtf="/mnt/tchandra-lab/Neil/genomes/GRCh38/Homo_sapiens.GRCh38.94.gtf"
export index="/mnt/tchandra-lab/Neil/genomes/${genome_build}/INDEXES/${genome_build}"
export threads=14

construct_sequences_folder.sh .

fastqc_module.sh . $threads

trimseqs_module.sh . $threads auto

hisat2_module.sh . $threads $index

processhisat2.sh . $threads $gtf

