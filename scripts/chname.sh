# !/bin/bash
# Cambiar el nombre a los ensambles
# bash chname.sh

cd /home/haydeeperuyero/MIEL/metagenomics/results/assembly/

for dname in *
do
	exname="${dname#assembly_}"
	mv ${dname}/final.contigs.fa ${dname}/${exname}_final.contigs.fa
done
