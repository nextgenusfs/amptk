
.. _quickstart:

AMPtk Quick Start
================

AMPtk can process NGS amplicon data from either 454, Illumina, or Ion Torrent sequencing platforms.
In order to get started, you will need to know how your amplicon reads were generated, i.e.
which primers you have used, PE or SE reads, where the amplicons have been barcoded, etc. For more
information see the :ref:`ReadLayout` information.

All AMPtk commands show a help menu when they are run without any parameters, for example::

    $ amptk

    Usage:       amptk <command> <arguments>
    version:     0.10.3

    Description: AMPtk is a package of scripts to process NGS amplicon data.  
                 Dependencies:  USEARCH v9.1.13 and VSEARCH v2.2.0
    
    Process:     ion         pre-process Ion Torrent data
                 illumina    pre-process folder of de-multiplexed Illumina data
                 illumina2   pre-process PE Illumina data from a single file
                 illumina3   pre-process PE Illumina + index reads (i.e. R1, R2, I)
                 454         pre-process Roche 454 (pyrosequencing) data
                 SRA         pre-process single FASTQ per sample data (i.e. SRA data)
             
    Clustering:  cluster     cluster OTUs (using UPARSE algorithm)
                 dada2       dada2 denoising algorithm (requires R, dada2, ShortRead)
                 unoise2     UNOISE2 denoising algorithm
                 cluster_ref closed/open reference based clustering (EXPERIMENTAL)

    Utilities:   filter      OTU table filtering
                 taxonomy    Assign taxonomy to OTUs
                 show        show number or reads per barcode from de-multiplexed data
                 select      select reads (samples) from de-multiplexed data
                 remove      remove reads (samples) from de-multiplexed data
                 sample      sub-sample (rarify) de-multiplexed reads per sample
                 drop        Drop OTUs from dataset
                 summarize   Summarize Taxonomy (create OTU-like tables and/or stacked bar graphs)
                 funguild    Run FUNGuild (annotate OTUs with ecological information) 
                 meta        pivot OTU table and append to meta data
                 heatmap     Create heatmap from OTU table
                 SRA-submit  De-multiplex data and create meta data for NCBI SRA submission

    Setup:       install     Download/install pre-formatted taxonomy DB. Only need to run once.
                 database    Format Reference Databases for Taxonomy
                 primers     List primers hard-coded in AMPtk. Can use in pre-processing steps.


