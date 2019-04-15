#!/bin/sh
echo "==>Arguments read:"
echo $@

#Running Kallisto inside quay.io/biocontainers/kallisto:0.45.0--hdcc98e5_0
cd /module
#echo "==>About to use Kallisto to create index"
# as of version 2.0 of this module, the index is part of the container
# kallisto index --index=data/index/Homo_sapiens.GRCh38.95_kalllisto_index data/index/Homo_sapiens.GRCh38.95.cdna.ncrna.VERSIONLESS.fa
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)."

echo "==>About to perform transcript quantitation."
# time kallisto quant --index=/build/index/Homo_sapiens.GRCh38.95_kalllisto_index --output-dir=RNASeq_quant --bias -b 2 data/fastqs/HT75CBCX2.1.ATTACTCG-ATAGAGGC/1.fastq.gz data/fastqs/HT75CBCX2.1.ATTACTCG-ATAGAGGC/2.fastq.gz data/fastqs/HT75CBCX2.2.ATTACTCG-ATAGAGGC/1.fastq.gz data/fastqs/HT75CBCX2.2.ATTACTCG-ATAGAGGC/2.fastq.gz
time kallisto quant --index=/build/index/Homo_sapiens.GRCh38.95_kalllisto_index --output-dir=RNASeq_quant $@

echo "==>About to perform transcript (TPM) aggregation."
time Rscript transcript2genes.R

echo "==>Done with all these analyses, the most inportant file created here is called 'gene_expression.csv'!"


# echo "==> DEBUGGING"
# echo "You may have dropped this from your clipboard:"
# echo "\t$pbpaste"
# echo "USE THIS URL:"
# echo "\t\thttp://127.0.0.1:8888/"
# jupyter-notebook --no-browser --port 8888 --ip=0.0.0.0 --NotebookApp.password='' --NotebookApp.token='' --NotebookApp.password_required=False --allow-root
