
.. _clustering:

AMPtk Clustering
================
Operational Taxonomic Unit (OTU) clustering has traditionally been done using a variety of methods, most simply clusters are generated based some sort of alignment between sequences and clusters are defined based on some percent identity, i.e. 97%.  More recently, the DADA2 denoising algorithm provides a clustering independent method that attempts to "correct" or "denoise" each sequence to a corrected sequence using statistical modeling of sequencing errors.  The clustering method you choose is highly dependent on what your research question is and what amplicons you are sequencing. 

All clustering scripts in AMPtk take the output of pre-processing (amptk ion, amptk illumina, etc) which is named ``output.demux.fg.gz``.  Once a run has been pre-processed, you can easily run several different clustering algorithms from the same input dataset. The first step of all AMPtk clustering scripts is to quality trim the reads using expected-errors quality trimming, which is literally an accumulation of the error probability for reach base pair in the read.  Thus, longer reads will have higher expected errors than shorter reads.  The quality trimmed reads are then pushed through each of the algorithms which output a set of OTUs.  In AMPtk, all OTU tables are generated from mapping the raw reads (no quality trimming done) to the OTUs (note you can change this parameter but you probably shouldn't). 

UPARSE clustering
-------------------------------------
The most commonly used method in OTU generation pipelines is UPARSE.  You can run UPARSE clustering in AMPtk as follows:

.. code-block:: none

    #run default settings
    amptk cluster -i mydata.demux.fq.gz -o mydata 
    
    #run with increased expected errors filter, perhaps you have really long reads
    amptk cluster -i mydata.demux.fq.gz -e 2.0 -o mydata

This script will generated a log file as well as OTUs named ``mydata.cluster.otus.fa`` and an OTU table ``mydata.otu_table.txt``.

DADA2 denoising
-------------------------------------
DADA2 can also be run from AMPtk, which utilizes a slightly modified version of DADA2 to retain the AMPtk data structure.

.. code-block:: none

    #run DADA2 on illumina dataset
    amptk dada2 -i mydata.demux.fq.gz --platform illumina -o mydata
    
    #run DADA2 on Ion Torrent dataset
    amptk dada2 -i mydata.demux.fq.gz --platform ion -o mydata
    
You will notice that the output of DADA2 results in 2 sets of output, the first is based on inferred sequences (iSeqs) or exact sequence variants (ESVs) and the second set out output comes from clustering the iSeqs into biological OTUs using UCLUST.  Depending on your research question, one or both of these outputs maybe useful.

UNOISE2 denoising
-------------------------------------
UNOISE2 is a denoising algorithm run in USEARCH9, which is similar in idea to DADA2. Likewise, the output data from UNOISE2 is the same as DADA2 -> where you get both multi-fasta and OTU tables for the iSeqs as well as clustered iSeqs.

.. code-block:: none

    #run UNOISE2
    amptk unoise2 -i mydata.demux.fq.gz -o mydata
    
UNOISE3 denoising
-------------------------------------
UNOISE3 is the denoising algorithm run in USEARCH10 (apparently the successor of UNOISE2). Likewise, the output data from UNOISE3 is the same as DADA2 and UNOISE2, although "it works better".... Note, for unoise3 you must have USEARCH10 installed.

.. code-block:: none

    #run UNOISE3
    amptk unoise3 -i mydata.demux.fq.gz -o mydata
