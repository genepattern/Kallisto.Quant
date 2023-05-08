#OLD #from jupyter/datascience-notebook:5ed91e8e3249
from jupyter/datascience-notebook:r-4.0.3
MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C
USER root
RUN apt-get update && apt-get install -y time

# Adding things to /build/
RUN mkdir /build

#Add Kallisto
# ADD BuildIndex/kallisto_linux-v0.46.1 /build/kallisto
# ADD BuildIndex/to_add /build/kallisto
ENV PATH="/opt/gtk/bin:${PATH}"

# Install sleuth: https://anaconda.org/bioconda/r-sleuth
RUN conda install -c bioconda r-sleuth==0.30.0
RUN conda install -c bioconda gtfparse==1.2.1


# Not used for now 2021-06-02 # Install tximport: https://anaconda.org/bioconda/bioconductor-tximport
#RUN conda install -c bioconda bioconductor-tximport==1.20.0
# Not used for now 2021-06-02 # Install biomaRt: 
#RUN conda install -c bioconda bioconductor-biomart==2.48.0

#Run this if you want to update indexes without changing the version of Kallisto, Sleuth,tximport, or biomaRt
#RUN wget -P /build/for_sleuth https://datasets.genepattern.org/data/module_support_files/Kallisto/for_sleuth/transcript2genes_dataframe.RDA
#RUN wget -P /build/index https://datasets.genepattern.org/data/module_support_files/Kallisto/index/Homo_sapiens.GRCh38.95_kalllisto_index
#RUN wget -P /build/index https://datasets.genepattern.org/data/module_support_files/Kallisto/index/Homo_sapiens.GRCh38.95.cdna.ncrna.VERSIONLESS.fa
#RUN wget -P /build/index https://datasets.genepattern.org/data/module_support_files/Kallisto/index/strip_ENSEMBL_versions.py

#Run this last
ADD src/ /module

# using the same user as the parent container as described in:
# https://hub.docker.com/r/jupyter/datascience-notebook/dockerfile
# USER $NB_UID

# build using this:
# docker build -t genepattern/kallisto:5.46.1 .
