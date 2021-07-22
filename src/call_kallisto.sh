#!/bin/bash

#Running Kallisto inside quay.io/biocontainers/kallisto:0.45.0--hdcc98e5_0

# for GP, you should not cd into /module, I think
# cd /module

#check if second argument ends with list.txt
echo $PATH
PATH="/build/kallisto:${PATH}"
echo $PATH

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

#Do somthing fancier here, if human, set INDEX to /build/kallisto/HUMAN_gencode.v37.transcripts.idx AND GENE_IDS to /build/kallisto/HUMAN_gencode.v37_t2g.csv
if [[ $1 == "Human" ]]
then
  INDEX=/build/kallisto/HUMAN_gencode.v37.transcripts.idx
  GENE_IDS=/build/kallisto/HUMAN_gencode.v37_t2g.csv
  echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS"
  echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS" >> module_log.txt
else
   if [ $1 == "Mouse" ]
   then
     INDEX=/build/kallisto/MOUSE_gencode.vM26.transcripts.idx
     GENE_IDS=/build/kallisto/MOUSE_gencode.vM26_t2g.csv
     echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS"
     echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS" >> module_log.txt
   else
     echo "We only have Human or Mouse indices at the moment"
     echo "We only have Human or Mouse indices at the moment" >> module_log.txt
   fi
fi
shift

OUT_BASENAME=$1
shift

#echo "==>About to use Kallisto to create index"
# as of version 2.0 of this module, the index is part of the container
# kallisto index --index=data/index/Homo_sapiens.GRCh38.95_kalllisto_index data/index/Homo_sapiens.GRCh38.95.cdna.ncrna.VERSIONLESS.fa
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)."
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)." >> module_log.txt

echo "==>About to perform transcript quantitation."
echo "==>About to perform transcript quantitation." >> module_log.txt
# NEW: # hardcoding "-b 2" [two bootstraps, the minimum, we don't use them in this workflow]
kallisto quant --index=$INDEX --output-dir=RNASeq_quant -b 2 $@ $list_of_files

echo "==>About to perform transcript (TPM) aggregation."
echo "==>About to perform transcript (TPM) aggregation." >> module_log.txt

Rscript /module/transcript2genes.R $GENE_IDS $OUT_BASENAME

mv RNASeq_quant/abundance.h5 "RNASeq_quant/${OUT_BASENAME}_abundance.h5"
mv RNASeq_quant/abundance.tsv "RNASeq_quant/${OUT_BASENAME}_abundance.tsv"

echo "==>Done with all these analyses, the most important file created here is called '$OUT_BASENAME.csv'!"
echo "==>Done with all these analyses, the most important file created here is called '$OUT_BASENAME.csv'!" >> module_log.txt


# echo "==> DEBUGGING"
# echo "You may have dropped this from your clipboard:"
# echo "\t$pbpaste"
# echo "USE THIS URL:"
# echo "\t\thttp://127.0.0.1:8888/"
# jupyter-notebook --no-browser --port 8888 --ip=0.0.0.0 --NotebookApp.password='' --NotebookApp.token='' --NotebookApp.password_required=False --allow-root
