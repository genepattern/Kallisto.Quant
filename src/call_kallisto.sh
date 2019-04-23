#!/bin/sh

#Running Kallisto inside quay.io/biocontainers/kallisto:0.45.0--hdcc98e5_0

# for GP, you should not cd into /module, I think
# cd /module

#check if second argument ends with list.txt

echo $PWD
echo "==>Arguments read:"
echo $@

echo $PWD > module_log.txt
echo "==>Arguments read:" >> module_log.txt
echo $@ >> module_log.txt

# Reading each line in the file list
filename=$1
list_of_files=""
while read line; do
  # reading each line
  #echo $line
  list_of_files="$list_of_files $line"
done < $filename
echo "==>list of files provided:"
echo "==>list of files provided:" >> module_log.txt
echo $list_of_files
echo $list_of_files >> module_log.txt

#==> parse $1 to open file, and turn newlines into spaces (create FILE_LIST_TO_STRING)
shift

#echo "==>About to use Kallisto to create index"
# as of version 2.0 of this module, the index is part of the container
# kallisto index --index=data/index/Homo_sapiens.GRCh38.95_kalllisto_index data/index/Homo_sapiens.GRCh38.95.cdna.ncrna.VERSIONLESS.fa
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)."
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)." >> module_log.txt

echo "==>About to perform transcript quantitation."
echo "==>About to perform transcript quantitation." >> module_log.txt
# OLD: # time kallisto quant --index=/build/index/Homo_sapiens.GRCh38.95_kalllisto_index --output-dir=RNASeq_quant --bias -b 2 data/fastqs/HT75CBCX2.1.ATTACTCG-ATAGAGGC/1.fastq.gz data/fastqs/HT75CBCX2.1.ATTACTCG-ATAGAGGC/2.fastq.gz data/fastqs/HT75CBCX2.2.ATTACTCG-ATAGAGGC/1.fastq.gz data/fastqs/HT75CBCX2.2.ATTACTCG-ATAGAGGC/2.fastq.gz
# NEW: # hardcoding "-b 2" [two bootstraps, the minimum, we don't use them in this workflow]
kallisto quant --index=/build/index/Homo_sapiens.GRCh38.95_kalllisto_index --output-dir=RNASeq_quant -b 2 $@ $list_of_files

echo "==>About to perform transcript (TPM) aggregation."
echo "==>About to perform transcript (TPM) aggregation." >> module_log.txt

Rscript /module/transcript2genes.R

echo "==>Done with all these analyses, the most inportant file created here is called 'gene_expression.csv'!"
echo "==>Done with all these analyses, the most inportant file created here is called 'gene_expression.csv'!" >> module_log.txt


# echo "==> DEBUGGING"
# echo "You may have dropped this from your clipboard:"
# echo "\t$pbpaste"
# echo "USE THIS URL:"
# echo "\t\thttp://127.0.0.1:8888/"
# jupyter-notebook --no-browser --port 8888 --ip=0.0.0.0 --NotebookApp.password='' --NotebookApp.token='' --NotebookApp.password_required=False --allow-root
