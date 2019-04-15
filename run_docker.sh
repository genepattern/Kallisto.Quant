# docker run -it --rm -v $PWD:/module -e GRANT_SUDO=yes genepattern/salmon:2.0 /module/Salmon_build_index.sh
# docker run -e TZ=America/Los_Angeles -it -p 8888:8888 --rm -v $PWD:/module -e GRANT_SUDO=yes genepattern/kallisto:1.0 /module/call_kallisto.sh

#from https://wurmlab.github.io/genomicscourse/2016-SIB/practicals/rnaseq/TP1
# For single-end data, the fragment length and standard deviation cannot be estimated directly from the data. The user needs to supply it (beware, fragment length is not read length, see https://groups.google.com/forum/#!topic/kallisto-sleuth-users/h5LeAlWS33w). This information has to be read from the Bioanalyzer/Fragment Analyzer results on the prepared RNA-seq libraries. For this practical, in the absence of this information, you will use length=200bp and sd=30, which should be close enough to real values.

docker run -e TZ=America/Los_Angeles -it --rm -v $PWD/module:/module -e GRANT_SUDO=yes genepattern/kallisto:3.0 time /module/call_kallisto.sh --single --fragment-length=200 --sd=30 -b 2 /module/test_data/SRR1515119_50k.fastq.gz
