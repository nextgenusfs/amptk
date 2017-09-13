.. AMPtk documentation master file, created by
   sphinx-quickstart on Tue Sep 12 09:56:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AMPtk documentation
=================================

.. toctree::
   :hidden:
   
   quick-start
   overview
   pre-processing
   clustering
   filtering
   taxonomy

AMPtk is a series of scripts to process NGS amplicon data using USEARCH and VSEARCH, it can also be used to process any NGS amplicon data and includes databases setup for analysis of fungal ITS, fungal LSU, bacterial 16S, and insect COI amplicons. It can handle Ion Torrent, MiSeq, and 454 data. At least USEARCH v9.1.13 and VSEARCH v2.2.0 are required as of AMPtk v0.7.0.

Install
==================
There are several ways to install AMPtk, you can download a `release <https://github.com/nextgenusfs/amptk/releases>`_. You can also build the latest unreleased version from github:
::
    git clone https://github.com/nextgenusfs/amptk.git
    export PATH="/path/to/amptk:$PATH"

But the easiest way to install AMPtk and its dependencies is with `HomeBrew <https://brew.sh>`_. For example:
::
    #install homebrew
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    
    #setup homebrew and link repositories
    brew doctor; brew tap homebrew/science; brew tap nextgenusfs/tap
    
    #install python dependencies
    pip install -U biopython natsort pandas numpy matplotlib edlib biom-format psutil
    
    #install AMPtk and dependencies
    brew install amptk
    
Dependencies
==================
1) AMPtk utilizes USEARCH9 which must be installed manually from the developer `here <http://www.drive5.com/usearch/download.html>`_.  Obtain the proper version of USEARCH v9.2.64 and softlink into the PATH:
::
    #make executable
    sudo chmod +x /path/to/usearch9.2.64_i86osx32
    
    #create softlink
    sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9

2) AMPtk requires VSEARCH, which you can install from `here <https://github.com/torognes/vsearch>`_. Note, if you use homebrew recipe it will be install automatically.
::
    #install vsearch with homebrew
    brew install vsearch
    
    #or with bioconda
    conda config --add channels bioconda
    conda install --yes vsearch

3) Several Python modules are also required, they can be installed with pip or conda:
::
    #install with pip
    pip install -U biopython natsort pandas numpy matplotlib edlib biom-format psutil
    
    #install with conda
    conda install biopython natsort pandas numpy matplotlib edlib biom-format psutil

4) DADA2 denoising algorithm requires installation of R and DADA2.  Instructions are located `here <http://benjjneb.github.io/dada2/>`_.
::
    #install with conda/bioconda
    conda config --add channels r
    conda install --yes bioconductor-dada2
    
More Information
==================

* :ref:`AMPtk Quick Start <quick-start>` - walkthrough of test data.
* :ref:`AMPtk Overview <overview>` - an overview of the steps in AMPtk.
* :ref:`AMPtk Pre-Processing <pre-processing>` - details of the critical pre-processing steps.
* :ref:`AMPtk Clustering <clustering>` - overview of clustering/denoising algorithms in AMPtk
* :ref:`AMPtk OTU Table Filtering <filtering>` - OTU table filtering based on Mock community
* :ref:`AMPtk Taxonomy <taxonomy>` - assigning taxonomy in AMPtk
