#!/bin/sh
# This program will perform taxonomic assignment with Kraken 2
# Use: it requires a Kraken database
# database: /data2/haydee/kraken_db

# Redirect output and error
exec 1>>  /data2/MIEL/results/taxonomy/kraken_PFP/output_pf.txt 2>> /data2/MIEL/results/taxonomy/kraken_PFP/error_pf.txt

# Moving to folder
if [ ! -d "/data2/MIEL/data/trimmed_fastq" ]; then
    echo "Error: El directorio /data2/MIEL/data/trimmed_fastq no existe." >&2
    exit 1
fi
cd /data2/MIEL/data/trimmed_fastq

# Assignation
kdat=$1
if [ -z "$kdat" ]; then
    echo "Error: No se proporcionó la base de datos de Kraken." >&2
    exit 1
fi

# Verify if files exist
ls *_1.trim.fastq.gz || { echo "No se encontraron archivos de entrada."; exit 1; }

# Taxonomic assignment with Kraken2
for file in *_1.trim.fastq.gz
do
    base=$(basename ${file} _1.trim.fastq.gz)
    echo "Procesando: ${base}" | tee -a output_pf.txt
    kraken2 --db $kdat --threads 12 --paired ${file} ${base}_2.trim.fastq.gz --output /data2/MIEL/results/taxonomy/kraken_PFP/krakens/${base}.kraken --report /data2/MIEL/results/taxonomy/kraken_PFP/reports/${base}.report
done
