#docker config for AMPtk installation 
FROM jupyter/scipy-notebook:latest
MAINTAINER Jon Palmer <nextgenusfs@gmail.com>

USER root

RUN apt-get update && apt-get install -y zlib1g-dev libssl-dev \
    libcurl4-gnutls-dev libxml2-dev libxslt1-dev

RUN conda config --add channels r && \
    conda config --add channels bioconda && \
    conda install --yes -p $CONDA_DIR/envs/python2 cython biopython natsort psutil \
    biom-format sra-tools rpy2 r-base r-essentials r-curl r-irkernel vsearch bioconductor-phyloseq

RUN pip2 install -U srapy

RUN ln -s /bin/tar /bin/gtar

RUN echo 'source("http://bioconductor.org/biocLite.R")' > ~/work/dada2_install.R && \
    echo 'biocLite(suppressUpdates = FALSE)' >> ~/work/dada2_install.R && \
    echo 'biocLite("ShortRead", suppressUpdates = FALSE)' >> ~/work/dada2_install.R && \
    echo 'biocLite("devtools")' >> ~/work/dada2_install.R && \
    echo 'library("devtools")' >> ~/work/dada2_install.R && \
    echo 'devtools::install_github("benjjneb/dada2")' >> ~/work/dada2_install.R

RUN $CONDA_DIR/envs/python2/bin/Rscript --vanilla ~/work/dada2_install.R

RUN rm ~/work/dada2_install.R

RUN git clone git://github.com/nextgenusfs/amptk.git && \
    cd amptk && \
    git checkout 45064ac3c54918410bc18a0a8c9fe3b234fc832d && \
    make && \
    cd ..

ENV PATH=/work:~/work/amptk:$PATH

RUN mkdir /work

WORKDIR /work