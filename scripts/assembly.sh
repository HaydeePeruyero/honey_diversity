# Para ensamblar los genomas
# Uso: bash assembly.sh
# Solo necesitamos estar 

# Crear folder

mkdir -p /files2/MIEL/results/assembly/

# Redirigir errores y outputs

#exec 2>> /files2/MIEL/results/assembly/output_assembly.txt 1>> /files2/MIEL/results/assembly/error_assembly.txt

# Moverse a la carpeta de datos

cd /files2/MIEL/data/trimmed_fastq2/

for file in *_1.trim.2.fastq.gz
do
	base=$(basename ${file} _1.trim.2.fastq.gz)
	echo "megahit -1 ${file} -2 ${base}_2.trim.2.fastq.gz -t 12 -o /files2/MIEL/results/assembly/assembly_${base}"
done
