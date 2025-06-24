# Crear los bioms por especie y por todos los archivos
# Uso: bash kraken_bioms.sh 

# Guardar error y output

exec 2>> /data2/MIEL/results/taxonomy/bioms/data_PFP/output_or.txt 1>> /data2/MIEL/results/taxonomy/bioms/data_PFP/error_or.txt

mkdir -p /data2/MIEL/results/taxonomy/bioms/data_PFP

cd /data2/MIEL/results/taxonomy/kraken_PFP/reports/

#kraken-biom *_Scapto.report --fmt json -o /files2/MIEL/results/taxonomy/bioms/db_original/Scapto.biom
#kraken-biom *_Melli.report --fmt json -o /files2/MIEL/results/taxonomy/bioms/db_original/Melli.biom
kraken-biom *.report --fmt json -o /data2/MIEL/results/taxonomy/bioms/data_PFP/honey_PFP.biom
