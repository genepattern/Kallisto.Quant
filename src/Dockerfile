from jupyter/datascience-notebook:5ed91e8e3249
MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C

#install Kallisto
RUN conda install -c bioconda kallisto==0.45.0

# Install sleuth
RUN conda install -c bioconda r-sleuth==0.30.0

# Install txtimport
RUN conda install -c bioconda bioconductor-tximport==1.8.0

# Install biomaRt
RUN conda install -c bioconda bioconductor-biomart==2.38.0
#RUN R -e 'install.packages("BiocManager")'
#RUN R -e 'BiocManager::install("biomaRt", version = "3.8")'

### Version 2 & 3

USER root
RUN apt-get update && apt-get install -y time

# old way
# ADD to_add /build/

# Adding things to /build/
RUN mkdir /build
RUN wget -P /build/for_sleuth https://datasets.genepattern.org/data/module_support_files/Kallisto/for_sleuth/transcript2genes_dataframe.RDA
RUN wget -P /build/index https://datasets.genepattern.org/data/module_support_files/Kallisto/index/Homo_sapiens.GRCh38.95_kalllisto_index
RUN wget -P /build/index https://datasets.genepattern.org/data/module_support_files/Kallisto/index/Homo_sapiens.GRCh38.95.cdna.ncrna.VERSIONLESS.fa
RUN wget -P /build/index https://datasets.genepattern.org/data/module_support_files/Kallisto/index/strip_ENSEMBL_versions.py

# RUN mkdir /module && cd /module
ADD . /module

# using the same user as the parent container as described in:
# https://hub.docker.com/r/jupyter/datascience-notebook/dockerfile
# USER $NB_UID

# build using this:
# docker build -t genepattern/kallisto:4.2 .
