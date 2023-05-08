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
echo "fastq_files: $1"
shift

#Do somthing fancier here, if human, set INDEX to /build/kallisto/HUMAN_gencode.v37.transcripts.idx AND GENE_IDS to /build/kallisto/HUMAN_gencode.v37_t2g.csv
echo "transcriptome: $1"
if [[ $1 == *"Homo_sapiens"* ]]
then
  INDEX=$1
  wget -P . $INDEX

  if [ -d ${INDEX##*/} ]
  then
    INDEX=${INDEX##*/}
  fi

  ls -alrt
  wget -P . https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
  python3 /module/run_t2g.py --transcript_gtf gencode.v43.annotation.gtf.gz
  GENE_IDS=HUMAN_gencode.v43_t2g.csv
  echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS"
  echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS" >> module_log.txt
elif [[ $1 == *"Mus_musculus"* ]]
then
  INDEX=$1
  wget -P . $INDEX
  
  if [ -d ${INDEX##*/} ]
  then
    INDEX=${INDEX##*/}
  fi

  wget -P . https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz
  python3 /module/run_t2g.py --transcript_gtf gencode.vM32.annotation.gtf.gz
  GENE_IDS=MOUSE_gencode.vM32_t2g.csv
  echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS"
  echo "Using INDEX=$INDEX and GENE_IDS=$GENE_IDS" >> module_log.txt
fi
shift

echo "bias: $1"
if [[ $1 == "YES" ]]
then
    BIAS="--bias"
else
    BIAS=""
fi
shift

echo "output_filename: $1"
OUT_BASENAME=$1
shift

echo "bootstrap_samples: $1"
BOOTSTRAP=$1
shift

echo "seed: $1"
SEED=$1
shift

echo "quantify_reads: $1"
if [[ $1 == "YES" ]]
then
    SINGLE="--single"
else
    SINGLE=""
fi
shift

echo "include_overhang: $1"
if [[ $1 == "YES" ]]
then
    OVERHANG="--single-overhang"
else
    OVERHANG=""
fi
shift

echo "read_direction: $1"
if [[ $1 == "forward" ]]
then
    DIRECTION="--fr-stranded"
elif [[ $1 == "reverse" ]]
then
    DIRECTION="--rf-stranded"
else
    DIRECTION=""
fi
shift

EXTRACTED_NUMBER=$(echo $1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' | tr '\n' ' '; echo "")
echo "fragment_length: $EXTRACTED_NUMBER"
if [[ ! -z "$EXTRACTED_NUMBER" ]]
then
    F_LENGTH="--fragment-length $EXTRACTED_NUMBER"
else
    F_LENGTH=""
fi
shift

EXTRACTED_NUMBER=$(echo $1 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' | tr '\n' ' '; echo "")
echo "fragment_sd: $EXTRACTED_NUMBER"
if [[ ! -z "$EXTRACTED_NUMBER" ]]
then
    F_SD="--sd $EXTRACTED_NUMBER"
else
    F_SD=""
fi
shift

echo "pseudobam: $1"
if [[ $1 == "YES" ]]
then
    PSEUDOBAM="--pseudobam"
else
    PSEUDOBAM=""
fi
shift

echo "genomebam: $1"
if [[ $1 == "YES" ]]
then
    GENOMEBAM="--genomebam"
else
    GENOMEBAM=""
fi
shift

EXTRACTED_PATH=$(echo "$1" | sed -e "s/^gtf_//")
echo "gtf_file: $EXTRACTED_PATH"
if [[ ! -z "$EXTRACTED_PATH" ]]
then
    GTF="--gtf $EXTRACTED_PATH"
else
    GTF=""
fi
shift

EXTRACTED_PATH=$(echo "$1" | sed -e "s/^chr_//")
echo "chromosome_file: $EXTRACTED_PATH"
if [[ ! -z "$EXTRACTED_PATH" ]]
then
    CHR="--chromosomes $EXTRACTED_PATH"
else
    CHR=""
fi
shift

#echo "==>About to use Kallisto to create index"
# as of version 2.0 of this module, the index is part of the container
# kallisto index --index=data/index/Homo_sapiens.GRCh38.95_kalllisto_index data/index/Homo_sapiens.GRCh38.95.cdna.ncrna.VERSIONLESS.fa
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)."
echo "==>Note that we are not calling Kallisto to create an index because for this application we have already built the index (it is part of the container)." >> module_log.txt

echo "==>About to perform transcript quantitation."
echo "==>About to perform transcript quantitation." >> module_log.txt
# NEW: # hardcoding "-b 2" [two bootstraps, the minimum, we don't use them in this workflow]
/module/kallisto/kallisto quant --index=$INDEX --output-dir=RNASeq_quant $BIAS --bootstrap-samples $BOOTSTRAP --seed $SEED $SINGLE $OVERHANG $DIRECTION $F_LENGTH $F_SD $PSEUDOBAM $GENOMEBAM $GTF $CHR $list_of_files

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
