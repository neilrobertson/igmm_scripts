#!/bin/sh

export PATH="/home/nerobert/bin/script_modules:$PATH"

#Set project specific params
export genome_build="GRCh38"
export gtf="/mnt/tchandra-lab/Neil/genomes/10x_CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
export threads=14

construct_sequences_folder.sh .

fastqc_module.sh . $threads

trimseqs_module.sh . $threads auto

star_module.sh . $threads $genome_build $gtf

STAR_summary_module.sh .
