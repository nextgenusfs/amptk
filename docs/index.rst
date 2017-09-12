.. AMPtk documentation master file, created by
   sphinx-quickstart on Tue Sep 12 09:56:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AMPtk documentation
=================================

.. toctree::
   :hidden:
   
   quick-start

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
    #install AMPtk
    brew install amptk

AMPtk utilizes USEARCH9 which must be installed manually from the developer `here <http://www.drive5.com/usearch/download.html>_.  Obtain the proper version of USEARCH v9.2.64 and softlink into the PATH:
::
    #make executable
    sudo chmod +x /path/to/usearch9.2.64_i86osx32
    #create softlink
    sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9

AMPtk also requires VSEARCH, which you can install from `here <>_. Note, if you use homebrew recipe it will be install automatically.

More Information
==================

* :ref:`quick-start`
* :ref:`overview`
* :ref:`pre-processing`
* :ref:`clustering`
* :ref:`filtering`
* :ref:`taxonomy`
