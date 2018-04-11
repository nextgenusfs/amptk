FROM continuumio/miniconda
MAINTAINER Jon Palmer <nextgenusfs@gmail.com>

RUN apt-get update && \
	apt-get install --fix-missing -y build-essential zlib1g-dev libssl-dev unzip && \
	export TAR=/bin/tar

RUN conda update conda && \
	conda clean --packages && \
	conda config --add channels r && \
	conda config --add channels defaults && \
	conda config --add channels conda-forge && \
	conda config --add channels bioconda && \
	conda install -y amptk r-tidyverse r-openssl

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
	Rscript -e "install.packages('devtools')" && \
	Rscript -e "options(unzip = 'internal'); devtools::install_github('tobiasgf/lulu')"

RUN mkdir /work
