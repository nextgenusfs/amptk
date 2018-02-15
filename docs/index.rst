.. AMPtk documentation master file, created by
   sphinx-quickstart on Tue Sep 12 09:56:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AMPtk documentation
=================================

.. toctree::
   :hidden:
  
   overview
   quick-start
   file-formats
   pre-processing
   clustering
   filtering
   taxonomy
   commands
   downstream


AMPtk is a series of scripts to process NGS amplicon data using USEARCH and VSEARCH, it can also be used to process any NGS amplicon data and includes databases setup for analysis of fungal ITS, fungal LSU, bacterial 16S, and insect COI amplicons. It can handle Ion Torrent, MiSeq, and 454 data. At least USEARCH v9.1.13 and VSEARCH v2.2.0 are required as of AMPtk v0.7.0.

Install
==================
There are several ways to install AMPtk, you can download a `release <https://github.com/nextgenusfs/amptk/releases>`_. You can also build the latest unreleased version from github:

.. code-block:: none

    git clone https://github.com/nextgenusfs/amptk.git
    export PATH="/path/to/amptk:$PATH"

But the easiest way to install AMPtk and its dependencies is with `HomeBrew <https://brew.sh>`_. For example:

.. code-block:: none

    #install homebrew
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    
    #setup homebrew and link repositories
    brew doctor; brew tap brewsci/bio; brew tap brewsci/science; brew tap nextgenusfs/tap
    
    #install python dependencies
    pip install -U biopython natsort pandas numpy matplotlib seaborn edlib biom-format psutil
    
    #install AMPtk and dependencies
    brew install amptk

You could also install all dependencies with Miniconda (conda). Download Miniconda2 for your operating system https://conda.io/miniconda.html

.. code-block:: none

    #install miniconda2 using installation script
    bash Miniconda2-latest-Linux-x86_64.sh
    
    #now setup your conda env with bioconda, type following in order
    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    
    #now install dependencies
    conda install vsearch biopython natsort pandas numpy matplotlib python-edlib \
        seaborn biom-format psutil r-base bioconductor-dada2 bioconductor-phyloseq
    
    #now get latest version of amptk
    git clone https://github.com/nextgenusfs/amptk.git
    export PATH="/path/to/amptk:$PATH"   

    
Dependencies
==================
1) AMPtk utilizes USEARCH9 which must be installed manually from the developer `here <http://www.drive5.com/usearch/download.html>`_.  Obtain the proper version of USEARCH v9.2.64 and softlink into the PATH:

.. code-block:: none

    #make executable
    sudo chmod +x /path/to/usearch9.2.64_i86osx32
    
    #create softlink
    sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9

1b) (optional) One script also requires USEARCH10, so you can download usearch10 and put into your path as follows:

.. code-block:: none

    #make executable
    sudo chmod +x /path/to/usearch10.0.240_i86osx32
    
    #create softlink
    sudo ln -s /path/to/usearch10.0.240_i86osx32 /usr/local/bin/usearch10

2) AMPtk requires VSEARCH, which you can install from `here <https://github.com/torognes/vsearch>`_. Note, if you use homebrew recipe it will be install automatically or can use conda.

.. code-block:: none

    #install vsearch with homebrew
    brew install vsearch
    
    #or with bioconda
    conda install -c bioconda vsearch

3) Several Python modules are also required, they can be installed with pip or conda:

.. code-block:: none

    #install with pip
    pip install -U biopython natsort pandas numpy matplotlib seaborn edlib biom-format psutil
    
    #install with conda
    conda install biopython natsort pandas numpy matplotlib seaborn python-edlib biom-format psutil

4) (optional)  DADA2 denoising algorithm requires installation of R and DADA2.  Instructions are located `here <http://benjjneb.github.io/dada2/>`_.

.. code-block:: none

    #install with conda/bioconda
    conda install r-base bioconductor-dada2

5) (optional) To run some preliminary community ecology stats via ``amptk stats`` you will also need the R package `Phyloseq <https://joey711.github.io/phyloseq/>`_.  One way to install with conda:

.. code-block:: none

    #install with conda/bioconda
    conda install r-base bioconductor-phyloseq
    
Run from Docker
==================
There is a base installation of AMPtk on Docker at nextgenusfs/amptk-base. Because usearch9 and usearch10 are required but must be personally licensed, here are the directions to get a working AMPtk docker image.

1) Download the Dockerfile build file.

.. code-block:: none

    wget https://raw.githubusercontent.com/nextgenusfs/amptk/master/Dockerfile

2) Download usearch9.2.64 and usearch10.0.240 for linux (32 bit version) `here <http://www.drive5.com/usearch/download.html>`_.


3) Build AMPtk docker image

.. code-block:: none

    docker build -t amptk -f Dockerfile .

4) You can now launch the docker image like so (make sure files you need are in current directory)

.. code-block:: none

    docker run -it --rm -v $PWD:/work amptk /bin/bash


More Information
==================

* :ref:`AMPtk Overview <overview>` - an overview of the steps in AMPtk.
* :ref:`AMPtk Quick Start <quick-start>` - walkthrough of test data.
* :ref:`AMPtk Pre-Processing <pre-processing>` - details of the critical pre-processing steps.
* :ref:`AMPtk Clustering <clustering>` - overview of clustering/denoising algorithms in AMPtk
* :ref:`AMPtk OTU Table Filtering <filtering>` - OTU table filtering based on Mock community
* :ref:`AMPtk Taxonomy <taxonomy>` - assigning taxonomy in AMPtk
* :ref:`AMPtk all commands <commands>` - all commands in AMPtk