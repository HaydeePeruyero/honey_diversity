# Cortar los archivos, remover los últimos 25 si no tienen calidad arriba de 28
for file in *_1.trim.fastq.gz
do
base=$(basename ${file} _1.trim.fastq.gz)
trimmomatic PE -threads 10 ${file} ${base}_2.trim.fastq.gz \
${base}_1.trim.2.fastq.gz ${base}_1un.trim.2.fastq.gz \
${base}_2.trim.2.fastq.gz ${base}_2un.trim.2.fastq.gz \
SLIDINGWINDOW:25:28 MINLEN:35 ILLUMINACLIP:TruSeq_LT_CD:2:40:15
done
