#!/bin/sh
# This program will perfom taxonomic assignation with Kraken 2
# Use: it need the kraken database
# Rafa database: /data/rafa/kraken-db
# Standard database /files/kraken/kraken-db
# nueva base PF /data2/haydee/kraken_db
# Redirect error and output

exec 2>> output_pf.txt 1>> error_pf.txt

# Moving to folder

cd /data2/MIEL/data/trimmed_fastq

#Assignations

kdat=$1

# Creation of folders

#mkdir -p /data2/MIEL/results/taxonomy/kraken_PF/krakens
#mkdir -p /data2/MIEL/results/taxonomy/kraken_PF/reports

# Taxonomic assignation with kraken2

for file in *_1.trim.fastq.gz
do
	base=$(basename ${file} _1.trim.fastq.gz)
	echo ${base}
	echo  kraken2 --db $kdat --threads 12 --paired ${file} ${base}_2.trim.fastq.gz --output /data2/MIEL/results/taxonomy/kraken_PF/krakens/${base}.kraken --report /data2/MIEL/results/taxonomy/kraken_PF/reports/${base}.report
done
