#!/bin/bash/
# Revisar la calidad de los bins encontrados
# Uso: bash check.sh
# Activar ambiente metagenomics y crear screen

cd /home/haydeeperuyero/MIEL/metagenomics/results/assembly/
for dir_name in *
do
extract_name="${dir_name#assembly_}"
mkdir -p ${dir_name}/CHECKM
checkm taxonomy_wf -t 25 domain Bacteria -x fasta ${dir_name}/MAXBIN/ ${dir_name}/CHECKM/
checkm qa ${dir_name}/CHECKM/Bacteria.ms ${dir_name}/CHECKM/ --file ${dir_name}/CHECKM/quality_${extract_name}.tsv --tab_table -o 2
done
