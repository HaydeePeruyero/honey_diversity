#!/bin/bash/
# Encontrar los Bins
# Uso: bash bins.sh
# Activar ambiente metagenomics y crear screen

cd /home/haydeeperuyero/MIEL/metagenomics/results/assembly/
for dir_name in *
do
extract_name="${dir_name#assembly_}"
mkdir -p ${dir_name}/MAXBIN
run_MaxBin.pl -thread 25 -contig ${dir_name}/${extract_name}_final.contigs.fa -reads /MIEL/data/trimmed_fastaq2/${extract_name}_1.trim.2.fastq.gz -reads2 /MIEL/data/trimmed_fastaq2/${extract_name}_2.trim.2.fastq.gz -out ${dir_name}/MAXBIN/${extract_name}
done
