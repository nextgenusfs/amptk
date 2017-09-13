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
   pre-processing
   clustering
   filtering
   taxonomy
   commands


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
    brew doctor; brew tap homebrew/science; brew tap nextgenusfs/tap
    
    #install python dependencies
    pip install -U biopython natsort pandas numpy matplotlib edlib biom-format psutil
    
    #install AMPtk and dependencies
    brew install amptk
    
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

2) AMPtk requires VSEARCH, which you can install from `here <https://github.com/torognes/vsearch>`_. Note, if you use homebrew recipe it will be install automatically.

.. code-block:: none

    #install vsearch with homebrew
    brew install vsearch
    
    #or with bioconda
    conda config --add channels bioconda
    conda install --yes vsearch

3) Several Python modules are also required, they can be installed with pip or conda:

.. code-block:: none

    #install with pip
    pip install -U biopython natsort pandas numpy matplotlib seaborn edlib biom-format psutil
    
    #install with conda
    conda install biopython natsort pandas numpy matplotlib seaborn edlib biom-format psutil

4) (optional)  DADA2 denoising algorithm requires installation of R and DADA2.  Instructions are located `here <http://benjjneb.github.io/dada2/>`_.

.. code-block:: none

    #install with conda/bioconda
    conda config --add channels r
    conda install --yes r-base bioconductor-dada2

5) (optional) To run some preliminary community ecology stats via ``amptk stats`` you will also need the R package `Phyloseq <https://joey711.github.io/phyloseq/>`_.  One way to install with conda:

.. code-block:: none

    #install with conda/bioconda
    conda config --add channels r
    conda install --yes r-base bioconductor-phyloseq
    
More Information
==================

* :ref:`AMPtk Overview <overview>` - an overview of the steps in AMPtk.
* :ref:`AMPtk Quick Start <quick-start>` - walkthrough of test data.
* :ref:`AMPtk Pre-Processing <pre-processing>` - details of the critical pre-processing steps.
* :ref:`AMPtk Clustering <clustering>` - overview of clustering/denoising algorithms in AMPtk
* :ref:`AMPtk OTU Table Filtering <filtering>` - OTU table filtering based on Mock community
* :ref:`AMPtk Taxonomy <taxonomy>` - assigning taxonomy in AMPtk
* :ref:`AMPtk all commands <commands>` - all commands in AMPtk