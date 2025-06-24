set -e # This will ensure that our script will exit if an error occurs
cd /files2/MIEL/data/

echo "Saving FastQC results..."
mv *.zip /files2/MIEL/results/fastqc_untrimmed_reads/
mv *.html /files2/MIEL/results/fastqc_untrimmed_reads/

cd /files2/MIEL/results/fastqc_untrimmed_reads/

echo "Unzipping..."
for filename in *.zip
    do
    unzip $filename
    done

echo "Saving summary..."
mkdir -p /files2/MIEL/docs
cat */summary.txt > /files2/MIEL/docs/fastqc_summaries.txt
