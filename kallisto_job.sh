#!/bin/sh

export PATH="/home/nerobert/bin/script_modules:$PATH"

#Set project specific params
export genome_build="GRCh38"
export gtf="/mnt/tchandra-lab/Neil/genomes/GRCh38/Homo_sapiens.GRCh38.94.chrfixed.gtf"
export index="/mnt/tchandra-lab/Neil/genomes/GRCh38/INDEXES/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx"
export threads=14

construct_sequences_folder.sh .

fastqc_module.sh . $threads

trimseqs_module.sh . $threads auto

conda activate
kallisto_module.sh . $threads $genome_build $index $gtf

Kallisto_summary_module.sh .
