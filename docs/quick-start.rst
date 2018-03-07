
.. _quick-start:

AMPtk Quick Start
================================

AMPtk can process NGS amplicon data from either 454, Illumina, or Ion Torrent sequencing platforms.
In order to get started, you will need to know how your amplicon reads were generated, i.e.
which primers you have used, PE or SE reads, where the amplicons have been barcoded, etc. For more
information see the :ref:`ReadLayout` information.

All AMPtk commands show a help menu when they are run without any parameters.

AMPtk comes with some very small test datasets, here we will go through a standard MiSeq Illumina test run. The data is located in the ``test_data`` folder of the amptk installation directory.  Move to that folder to begin.

Pre-Processing Reads
-------------------------------------
The first step is to pre-process or demultiplex your reads.  In this example, the data is PE MiSeq (2 x 300 bp) and it has been processed by the sequencing center using the Illumina program ``bcl2fastq`` so each sample has a forward and reverse read file (denoted by _R1 and _R2 in the sample name). The folder ``illumina_test_data`` contains 3 PE samples we will use for a demo.  Since the data is in this format, we will first use the ``amptk illumina`` format.  Note that some sequencing centers distribute their data in different formats, ``amptk illumina2`` and ``amptk illumina3`` are additional methods that can be used to pre-process illumina data.

.. code-block:: none

    #run the pre-processing script, specifying primers used
    amptk illumina -i illumina_test_data -o miseq -f fITS7 -r ITS4
    
    -------------------------------------------------------
    [03:53:16 PM]: OS: MacOSX 10.12.6, 8 cores, ~ 16 GB RAM. Python: 2.7.12
    [03:53:16 PM]: AMPtk v1.0.0, USEARCH v9.2.64, VSEARCH v2.4.4
    [03:53:16 PM]: Gzipped files detected, uncompressing
    [03:53:16 PM]: Merging Overlaping Pairs using USEARCH
    [03:53:16 PM]: working on sample 301-1 (Read Length: 300)
    [03:53:17 PM]: 100 reads passed (100.0%)
    [03:53:17 PM]: working on sample 301-2 (Read Length: 300)
    [03:53:17 PM]: 100 reads passed (100.0%)
    [03:53:17 PM]: working on sample spike (Read Length: 300)
    [03:53:17 PM]: 400 reads passed (100.0%)
    [03:53:17 PM]: Stripping primers and trim to 300 bp
    [03:53:17 PM]: splitting the job over 8 cpus, but this may still take awhile
    [03:53:17 PM]: Foward primer: GTGARTCATCGAATCTTTG,  Rev comp'd rev primer: GCATATCAATAAGCGGAGGA
    -------------------------------------------------------
    [03:53:18 PM]: Concatenating Demuxed Files
    [03:53:18 PM]: 600 total reads
    [03:53:18 PM]: 600 Fwd Primer found, 426 Rev Primer found
    [03:53:18 PM]: 0 discarded too short (< 100 bp)
    [03:53:18 PM]: 600 valid output reads
    [03:53:18 PM]: Found 3 barcoded samples
                            Sample:  Count
                             spike:  400
                             301-1:  100
                             301-2:  100
    [03:53:18 PM]: Output file:  miseq.demux.fq.gz (53.1 KB)
    [03:53:18 PM]: Mapping file: miseq.mapping_file.txt
    -------------------------------------------------------


You should see that script found 600 valid reads, and output them into a file called ``miseq.demux.fq.gz`` and you'll note that it also output a QIIME-like mapping file ``miseq.mapping_file.txt`` which can be used later on to add meta data to and generate a BIOM file (during assigning taxonomy step).  You can find a detailed log file in ``miseq.amptk-demux.log``.  The script has merged the PE reads using usearch, removed phiX spike-in, removed the forward and reverse primers, relabeled the FASTQ headers, and then concatenated all samples together.  Note that the default settings require that reads have a valid forward primer, if you used a custom sequencing primer then you will need to turn this setting off with ``--require_primer off``.

Clustering Data
-------------------------------------
We can now take the output here and run the UPARSE clustering algorithm.  The steps of UPARSE are to quality filter with expected errors, dereplicate, sort by abundance, cluster using 97% identity to create OTUs, map original reads back to OTUs to make OTU table.

.. code-block:: none

    amptk cluster -i miseq.demux.fq.gz -o miseq
    
    -------------------------------------------------------
    [03:54:29 PM]: OS: MacOSX 10.12.6, 8 cores, ~ 16 GB RAM. Python: 2.7.12
    [03:54:29 PM]: AMPtk v1.0.0, USEARCH v9.2.64, VSEARCH v2.4.4
    [03:54:29 PM]: Loading FASTQ Records
    [03:54:29 PM]: 600 reads (53.1 KB)
    [03:54:29 PM]: Quality Filtering, expected errors < 1.0
    [03:54:29 PM]: 394 reads passed
    [03:54:29 PM]: De-replication (remove duplicate reads)
    [03:54:29 PM]: 112 reads passed
    [03:54:29 PM]: Clustering OTUs (UPARSE)
    [03:54:29 PM]: 32 OTUs
    [03:54:29 PM]: Cleaning up padding from OTUs
    [03:54:29 PM]: Mapping Reads to OTUs and Building OTU table
    [03:54:30 PM]: 429 reads mapped to OTUs (72%)
    -------------------------------------------------------
    OTU Clustering Script has Finished Successfully
    -------------------------------------------------------
    Clustered OTUs: /Users/jon/amptk/test_data/miseq.cluster.otus.fa
    OTU Table: /Users/jon/amptk/test_data/miseq.otu_table.txt
    -------------------------------------------------------

Filtering Data
-------------------------------------
Since we included a mock community in our sample, we will filter for index-bleed. Note in this toy example there is no index-bleed.

.. code-block:: none

    amptk filter -i miseq.otu_table.txt -f miseq.cluster.otus.fa -b spike -m mock2
    
    -------------------------------------------------------
    [03:56:53 PM]: OS: MacOSX 10.12.6, 8 cores, ~ 16 GB RAM. Python: 2.7.12
    [03:56:54 PM]: AMPtk v1.0.0, USEARCH v9.2.64, VSEARCH v2.4.4
    [03:56:54 PM]: Loading OTU table: miseq.otu_table.txt
    [03:56:54 PM]: OTU table contains 32 OTUs and 429 read counts
    [03:56:54 PM]: Mapping OTUs to Mock Community (USEARCH)
    [03:56:54 PM]: Sorting OTU table naturally
    [03:56:54 PM]: Removing OTUs according to --min_reads_otu: (OTUs with less than 2 reads from all samples)
    [03:56:54 PM]: Normalizing OTU table to number of reads per sample
    [03:56:54 PM]: Index bleed, mock into samples: 0.000000%.  Index bleed, samples into mock: 0.000000%.
    [03:56:54 PM]: No spike-in mock (-b) or index-bleed (-p) specified, thus not running index-bleed filtering
    [03:56:54 PM]: spike sample has 24 OTUS out of 26 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.040%
    [03:56:54 PM]: Filtering OTU table down to 8 OTUs and 136 read counts
    [03:56:54 PM]: Filtering valid OTUs
    -------------------------------------------------------
    OTU Table filtering finished
    -------------------------------------------------------
    OTU Table Stats:      miseq.stats.txt
    Sorted OTU table:     miseq.sorted.txt
    Normalized/filter:    miseq.normalized.txt
    Final Binary table:   miseq.final.binary.txt
    Final OTU table:      miseq.final.txt
    Filtered OTUs:        miseq.filtered.otus.fa
    -------------------------------------------------------

Post clustering LULU
-------------------------------------
AMPtk can run the LULU post clustering tool that can identify potential errors in a dataset. This step is optional.

.. code-block:: none

    amptk lulu -i miseq.final.txt -f miseq.filtered.otus.fa -o miseq

    -------------------------------------------------------
    [Mar 07 11:47 AM]: OS: MacOSX 10.13.3, 8 cores, ~ 16 GB RAM. Python: 2.7.14
    [Mar 07 11:47 AM]: AMPtk v1.1.0, USEARCH v9.2.64, VSEARCH v2.6.2
    [Mar 07 11:47 AM]: R v3.4.1; LULU v0.1.0
    [Mar 07 11:47 AM]: Loading 8 OTUs
    [Mar 07 11:47 AM]: Generating pairwise percent identity between OTUs using VSEARCH at 84% identity
    [Mar 07 11:47 AM]: Running LULU algorithm
    [Mar 07 11:47 AM]: LULU has merged 1 OTUs, output data contains 7 OTUs
    [Mar 07 11:47 AM]: LULU OTU table post processing finished
    ----------------------------------
    OTU table:  miseq.lulu.otu_table.txt
    OTU FASTA:  miseq.lulu.otus.fa
    LULU map:   miseq.lulu.otu-map.txt
    ----------------------------------


Assign Taxonomy
-------------------------------------
We can now assign taxonomy to our OTUs and create the final BIOM output file.

.. code-block:: none

    amptk taxonomy -f miseq.lulu.otus.fa -i miseq.lulu.otu_table.txt -m miseq.mapping_file.txt -d ITS2 -o miseq
    
    -------------------------------------------------------
    [Mar 07 12:51 PM]: OS: MacOSX 10.13.3, 8 cores, ~ 16 GB RAM. Python: 2.7.14
    [Mar 07 12:51 PM]: AMPtk v1.1.0, USEARCH v9.2.64, VSEARCH v2.6.2
    [Mar 07 12:51 PM]: Loading FASTA Records
    [Mar 07 12:51 PM]: 7 OTUs
    [Mar 07 12:51 PM]: Global alignment OTUs with usearch_global (USEARCH)
    [Mar 07 12:51 PM]: Classifying OTUs with UTAX (USEARCH)
    [Mar 07 12:51 PM]: Classifying OTUs with SINTAX (USEARCH)
    [Mar 07 12:52 PM]: Appending taxonomy to OTU table and OTUs
    [Mar 07 12:52 PM]: Generating phylogenetic tree
    [Mar 07 12:52 PM]: Taxonomy finished: miseq.taxonomy.txt
    [Mar 07 12:52 PM]: Classic OTU table with taxonomy: miseq.otu_table.taxonomy.txt
    [Mar 07 12:52 PM]: BIOM OTU table created: miseq.biom
    [Mar 07 12:52 PM]: OTUs with taxonomy: miseq.otus.taxonomy.fa
    [Mar 07 12:52 PM]: OTU phylogeny: miseq.tree.phy
    -------------------------------------------------------


And we now have an OTU table complete with taxonomy:

.. code-block:: none

    #OTU ID	301-1	301-2	spike	Taxonomy
    OTU1	0	57	0	GSL|100.0|KX857803;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Hymenochaetales,f:Schizoporaceae
    OTU10	0	0	16	GSL|100.0|KF169595;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Fomitopsidaceae,g:Fomitopsis
    OTU11	0	0	14	GS|100.0|KP055600;k:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,f:Vibrisseaceae,g:Phialocephala,s:Phialocephala lagerbergii
    OTU12	0	0	14	GS|100.0|KU668958;k:Fungi,p:Ascomycota
    OTU13	0	0	10	GS|99.7|EU402583;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Polyporaceae,g:Leptoporus,s:Leptoporus mollis
    OTU14	0	0	10	GSL|100.0|KP135060;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales
    OTU15	0	0	10	GS|100.0|KU668956;k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Eurotiales,f:Trichocomaceae,g:Penicillium,s:Penicillium nothofagi
    OTU16	0	0	13	GS|100.0|KM493837;k:Fungi
    OTU17	0	0	10	GS|100.0|KU668960;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Fomitopsidaceae,g:Laetiporus,s:Laetiporus caribensis
    OTU18	0	0	11	GS|100.0|KY886708;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Fomitopsidaceae,g:Laetiporus,s:Laetiporus cremeiporus
    OTU19	0	0	10	GS|100.0|KU668959;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Polyporaceae,g:Wolfiporia,s:Wolfiporia dilatohypha
    OTU2	36	0	0	GSL|100.0|HQ222028;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Strophariaceae,g:Pholiota
    OTU20	7	3	0	GS|100.0|AY465463;k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Chaetothyriales,f:Herpotrichiellaceae,g:Phialophora
    OTU21	0	0	10	GS|100.0|KU668966;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Phanerochaetaceae,g:Antrodiella,s:Antrodiella semisupina
    OTU22	0	0	9	GSL|100.0|AY340039;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Hymenochaetales,f:Hymenochaetaceae,g:Phellinus,s:Phellinus cinereus
    OTU23	0	0	9	GS|100.0|KM373239;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Polyporaceae,g:Trametes,s:Trametes gibbosa
    OTU24	0	0	6	GS|100.0|EU402560;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Fomitopsidaceae,g:Laetiporus,s:Laetiporus sulphureus
    OTU25	0	0	8	GS|100.0|KU668954;k:Fungi,p:Mortierellomycota
    OTU26	0	0	9	GS|100.0|DQ398958;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Corticiales,f:Punctulariaceae,g:Punctularia,s:Punctularia strigosozonata
    OTU27	19	0	0	GS|100.0|DQ647503;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Hymenochaetales,g:Peniophorella,s:Peniophorella pubera
    OTU28	0	0	10	GS|100.0|KU668961;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Fomitopsidaceae,g:Laetiporus,s:Laetiporus persicinus
    OTU29	0	4	0	US|0.8150|LN827700;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales
    OTU3	0	0	21	GS|100.0|KU668970;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Meruliaceae,g:Bjerkandera,s:Bjerkandera adusta
    OTU30	0	5	0	GSL|100.0|KF928438;k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Chaetothyriales,f:Herpotrichiellaceae,g:Exophiala,s:Exophiala cancerae
    OTU31	2	0	0	GS|100.0|KX221389;k:Fungi
    OTU32	1	2	0	GS|100.0|JX946684;k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Venturiales,f:Venturiaceae,g:Protoventuria,s:Protoventuria alpina
    OTU4	0	0	14	GS|100.0|KU668955;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Hymenochaetales,f:Schizoporaceae,g:Schizopora
    OTU5	0	0	15	GS|100.0|KU668973;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Phanerochaetaceae,g:Phanerochaete,s:Phanerochaete laevis
    OTU6	0	0	17	GS|100.0|KU668975;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Polyporaceae,g:Leptoporus,s:Leptoporus mollis
    OTU7	0	0	17	GS|100.0|KU668968;k:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales
    OTU8	0	0	17	GSL|100.0|KP216946;k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Eurotiales,f:Trichocomaceae,g:Aspergillus
    OTU9	0	0	13	GS|100.0|EU402553;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Polyporales,f:Fomitopsidaceae,g:Laetiporus,s:Laetiporus gilbertsonii

We then also have a BIOM file complete with any metadata and taxonomy (``miseq.biom``). This file could be taken to some of the :ref:`downstream processing software <downstream>`. 


