#!/usr/bin/env bash

export dir=$1
cd $dir

"Working in: ${dir}"

export PYTHONPATH="/mnt/tchandra-lab/Neil/repository/workspace/Baselib:$PYTHONPATH"

#Set project specific params
export genome_build="GRCh38"
export threads=10

export chrm_sizes="/mnt/tchandra-lab/Neil/genomes/${genome_build}/chrmSizes.${genome_build}"


mkdir peaks_SICER
cd peaks_SICER

IFS=$'\n' GLOBIGNORE='*' command eval  'MANIFEST=($(cat ../ChIP_comparisons.txt))'
for line in "${MANIFEST[@]}"; 
do
 IFS=$'\t'
 tmp=($line)
 export treatment_file="${tmp[0]}"
 export input_control="${tmp[1]}"
 export output_name="${tmp[2]}"

 echo "Working on: ${output_name} \n"
 mkdir ./peaks_SICER/${output_name}

 epic2  --treatment ../bam/${treatment_file}/${treatment_file}.uniq.rmdup.bed \
--control ../bam/${input_control}/${input_control}.uniq.rmdup.bed \
--genome hg38 \
--fragment-size 200 \
--effective-genome-fraction 0.75 \
--chromsizes ${chrm_sizes} \
-fdr 0.01 \
--keep-duplicates

done

rm -rf ./bed/

echo "Complete!!"