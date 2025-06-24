#!/bin/bash/
# Correr gtdbk
# Uso: activar ambiente gtbdk
# bash gtdbk_identify.sh
cd /files2/MIEL/results/assembly/

for dir_name in *
do
extract_name="${dir_name#assembly_}"
echo "mkdir -p /files2/MIEL/results/gtbd/${extract_name}/genomes"
echo "mkdir -p /files2/MIEL/results/gtdb/${extract_name}/align"
echo "mkdir -p /files2/MIEL/results/gtdb/${extract_name}/classify"
echo "mkdir -p /files2/MIEL/results/gtdb/${extract_name}/identify"
echo "cp ${dir_name}/MAXBIN/*.fasta gtdb/${extract_name}/genomes/."
echo "gtdbtk identify --genome_dir /files2/MIEL/results/gtdb/${extract_name}/genomes --out_dir /files2/MIEL/results/gtbd/${extract_name}/identify --extension fasta --cpus 2"
echo "gtdbtk align --identify_dir /files2/MIEL/results/gtdb/${extract_name}/identify --out_dir /files2/MIEL/results/gtdb/${extract_name}/align --cpus 2"
echo "gtdbtk classify --genome_dir /files2/MIEL/results/gtdb/${extract_name}/genomes --align_dir /files2/MIEL/results/gtdb/${extract_name}/align --out_dir /files2/MIEL/results/gtdb/${extract_name}/classify -x fasta --cpus 2"
done

