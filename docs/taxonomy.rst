
.. _taxonomy:

AMPtk Taxonomy
================

AMPtk can assign taxonomy using Blast, RDP Classifier, Global Alignment, UTAX, and/or SINTAX. The default method in AMPtk is 'hybrid' taxonomy assignment which calculates a consensus LCA (last common ancestor) taxonomy based on the results of Global Alignment (USEARCH/VSEARCH), UTAX, and SINTAX. This method is a conservative estimate of taxonomy because of the use of a LCA algorithm, however, is more robust than using any of the searches individually. The LCA method used here is applied twice, first in finding the LCA of identical Global Alignment hits and then again to compare the taxonomy from Global Alignment to the Bayesian Classifiers. If there is no conflict between taxonomies from each method, the taxonomy string with most information is retained.

**The hybrid taxonomy algorithm works as follows:**

1) Global Alignment is performed on a reference database. The method captures all of the top hits from the database that have Percent Identity higher than ``--usearch_cutoff`` (Default: 0.7). The top hits that have the most taxonomy information (levels) are retained and LCA is run on these hits.
2) UTAX Classifier is run generating a taxonomy string that scores better than ``--utax_cutoff`` (Default: 0.8).
3) SINTAX Classifier is run generating a taxonomy string that scores better than ``--sintax_cutoff`` (Default: 0.8).
4) The Bayesian Classifier results are compared, and the method that produces the most levels of taxonomy above the threshold is retained.
5) If the best Global Alignment result is less than 97% identical, then final taxonomy string defaults to the best Bayesian Classifier result (UTAX or SINTAX).
6) If the best Global Alignment result is greater than 97% identical then that hit is retained. A final LCA algorithm is applied to the Global Alignment hit and the Best Bayesian Classifier hit.

**The resulting taxonomy strings have the following format:**

.. code-block:: none

    #OTUID	taxonomy
    OTU1	GSL|100.0|KP776995;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Coniochaetales,f:Coniochaetaceae,g:Coniochaeta
    OTU2	US|0.9389|HF674734;k:Fungi,p:Ascomycota
    OTU3	GS|100.0|AY465463;k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Chaetothyriales,f:Herpotrichiellaceae,g:Phialophora
    OTU4	GS|98.0|KX222656;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales,f:Lasiosphaeriaceae
    OTU5	US|1.0000|FM200433;k:Fungi
    OTU6	GS|100.0|KF850373;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales,f:Chaetomiaceae,g:Trichocladium,s:Trichocladium opacum
    OTU7	GDL|100.0|KF660573;k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Pleosporales
    OTU8	GS|100.0|KX611531;k:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,f:Helotiaceae,g:Varicosporium,s:Varicosporium elodeae
    OTU9	GSL|100.0|JN655624;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Hypocreales,g:Ilyonectria
    OTU10	SS|0.9000|KT197416;k:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales

Taxonomy information follows the ; in the string. The pre-taxonomy identifiers are coded to mean the following:

.. code-block:: none

    #OTUID	taxonomy
    OTU1	Method|Pident/Score|Ref_DB ID;

Coding for the different methods:

- **GS**:    Global Alignment (G), Bayesian and Global Alignment Taxonomy is the Same (S)
- **GSL**:   Global Alignment (G), Bayesian and Global Alignment Taxonomy is the Same (S), LCA algorithm used (L)
- **GDL**:   Global Alignment (G), Bayesian and Global Alignment Taxonomy Discrepancy (D), LCA algorithm used (L)
- **US**:    UTAX Classifier (U), Bayesian and Global Alignment Taxonomy is the Same (S)
- **SS**:    SINTAX Classifier (S), Bayesian and Global Alignment Taxonomy is the Same (S)

Assigning Taxonomy
-------------------------------------
A typically command to assign taxonomy in AMPtk looks like this:

.. code-block:: none

    amptk taxonomy -i input.otu_table.txt -f input.cluster.otus.fa -m input.mapping_file.txt -d ITS2
    
This command will run the default hybrid method and will use the ITS2 database (``-d ITS2``).  The output of this command will be a tab delimited taxonomy file ``input.taxonomy.txt``, an OTU table contain taxonomy as the last column ``input.otu_table.taxonomy.txt```, a multi-fasta file with taxonomy as OTU headers ``input.otus.taxonomy.fa``, and BIOM file containing taxonomy and metadata ``input.biom``.
     
Taxonomy Databases
-------------------------------------
AMPtk is packaged with 4 reference databases for fungal ITS, fungal LSU, bacterial 16S, and arthropod/chordate mtCOI. These pre-built databases are updated frequently when reference databases are updated and can be downloaded/installed as follows:

.. code-block:: none

    #install all databases
    amptk install -i ITS 16S LSU COI

    #install only ITS database
    amptk install -i ITS

    #update database
    amptk install -i ITS 16S LSU COI --force
    
Users can also build their own custom databases, with the largest obstacle to overcome being formatting the taxonomy headers for reference databases.  Because AMPtk uses UTAX/SINTAX Bayesian classifiers, it uses the same taxonomy header formatting which looks like the following ``Kingdom(k), Phylum(p), Class(c), Order(o), Family(f), Genus(g), Species(s)``:

.. code-block:: none

    >BOLD:ACI6695;tax=k:Animalia,p:Arthropoda,c:Insecta,o:Coleoptera,f:Elateridae,g:Nipponoelater,s:Nipponoelater babai
    >S004604051;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Hymenochaetales,f:Hymenochaetaceae,g:Inonotus,s:Sanghuangporus zonatus
    >S004127186;tax=k:Fungi,p:Ascomycota
    >S004061552;tax=k:Fungi,p:Ascomycota,c:Eurotiomycetes,s:Pyrenula sanguinea

Note that if levels of taxonomy are unknown they can be left out, but should not contain things like `unclassified`, `unknown`, `incertae_sedis` -> as these levels of taxonomy are not informative and will produce undesired results.

Taxonomy databases are built with the ``amptk database`` command.  This command contains some parsers for known fasta header datasets, however, it is likely that generating custom databases will require some scripting to reformat the fasta headers.  The pre-build databases for AMPtk were constructed as follows:

**Fungal ITS DB**

These databases were created from Unite v7.2.2 (released June 28th, 2017), first downloading two databases from the UNITE website.  First the General FASTA release of the DB `here <https://unite.ut.ee/sh_files/sh_general_release_28.06.2017.zip>`_, and `here <https://unite.ut.ee/sh_files/sh_general_release_s_28.06.2017.zip>`_.  Then the Full UNITE+INSD database `here <https://unite.ut.ee/sh_files/UNITE_public_28.06.2017.fasta.zip>`_.  For the general FASTA releases, the 'developer' fasta files are used. The taxonomy information is then reformated and databases produced as follows:

.. code-block:: none

    #Create full length ITS USEARCH Database, convert taxonomy, and create USEARCH database
    amptk database -i UNITE_public_28.06.2017.fasta -o ITS -f ITS1-F -r ITS4 \
        --create_db usearch --keep_all

    #Create UTAX Databases
    amptk database -i sh_general_release_dynamic_28.06.2017_dev.fasta  \
        -o ITS_UTAX --create_db utax -f ITS1-F -r ITS4 --keep_all
        --derep_fulllength --lca
        
    amptk database -i sh_general_release_dynamic_s_28.06.2017_dev.fasta \
        -o ITS1_UTAX --create_db utax -f ITS1-F -r ITS2 --keep_all
        --derep_fulllength --lca
        
    amptk database -i sh_general_release_dynamic_s_28.06.2017_dev.fasta \
        -o ITS2_UTAX --create_db utax -f fITS7 -r ITS4 --derep_fulllength --lca

**Arthropod/Chordate mtCOI DB**

These data were pulled from the `BOLDv4 database <http://v4.boldsystems.org>`_  Since most studies using mtCOI regions are interested in identification of insects in the diets of animals, the BOLD database was queried as follows.  All Chordata sequences were downloaded by querying the `BIN database using the search term Chordata <http://v4.boldsystems.org/index.php/Public_BINSearch?query=Chordata&searchBIN=Search+BINs>`_.  Similarly, the Arthropods were searched by querying the `BIN databases using the search term Arthropoda <http://v4.boldsystems.org/index.php/Public_BINSearch?query=Arthropoda&searchBIN=Search+BINs>`_.  All data was then downloaded as TSV output.

The TSV output files (~ 6GB) where then each formatted using the following method, which reformats the taxonomy information and pulls sequences that are annotated in BINS and then clusters sequences in each bin to 99%.

.. code-block:: none

    #reformat taxonomy
    amptk/util/bold2utax.py -i Arthropoda_bold_data.txt -o arthropoda.bold.bins.fa
    amptk/util/bold2utax.py -i Chordata_bold_data.txt -o chordata.bold.bins.fa

    #combine datasets
    cat arthropoda.bold.bins.fa chordata.bold.bins.fa > all.data.bins.fa
    
    #generate global alignment database
    amptk database -i all.data.bins.fa --skip_trimming --keep_all --min_len 125 \
        --derep_fulllength --create_db usearch -o COI --format off

The data is then further processed with a second script that will search for priming sites and then randomly subsample the data down to a number of records that can be used to train UTAX and then database was created.

.. code-block:: none

 #searches for priming sites and subsamples to 90,000 records
 amptk/util/bold2amptk.py -i all.data.bins.fa -o arthropods.chordates
 
 #generate utax database
 amptk database -i arthropods.chordates.genus4utax.fa -o COI_UTAX \
    --format off --create_db utax --skip_trimming

**LSU database**

The fungal 28S database (LSU) was downloaded from `RDP <http://rdp.cme.msu.edu/download/current_Fungi_unaligned.fa.gz>`_.  The sequences were then converted into AMPtk databases as follows:

.. code-block:: none

 amptk database -i fungi.unaligned.fa -o LSU --format rdp2utax \
    --skip_trimming --create_db usearch --derep_fulllength --keep_all

To generate a training set for UTAX, the sequences were first dereplicated, and clustered at 97% to get representative sequences for training.  This training set was then converted to a UTAX database:

.. code-block:: none

 amptk database -i fungi.trimmed.fa -o LSU_UTAX --format off \
    --skip_trimming --create_db utax --keep_all

**16S database**
This is downloaded from `R. Edgar's website <http://drive5.com/utax/data/rdp_v16.tar.gz>`_ and then formatted for AMPtk.  Note there is room for substantial improvement here, I just don't typically work on 16S - so please let me know if you want some suggestions on what to do here.

.. code-block:: none

 amptk database -i rdp_v16.fa -o 16S --format off --create_db utax \
    --skip_trimming --keep_all

Checking Installed Databases
-------------------------------------
A simple ``amptk taxonomy`` command will show you all the arguments as well as display which databases have been installed.

.. code-block:: none

    amptk taxonomy

    Usage:       amptk taxonomy <arguments>
    version:     1.0.0

    Description: Script maps OTUs to taxonomy information and can append to an OTU table (optional).  
                 By default the script uses a hybrid approach, e.g. gets taxonomy information from 
                 SINTAX, UTAX, and global alignment hits from the larger UNITE-INSD database, and 
                 then parses results to extract the most taxonomy information that it can at 'trustable' 
                 levels. SINTAX/UTAX results are used if BLAST-like search pct identity is less than 97%.  
                 If % identity is greater than 97%, the result with most taxonomy levels is retained.
    
    Arguments:   -f, --fasta         Input FASTA file (i.e. OTUs from amptk cluster) (Required)
                 -i, --otu_table     Input OTU table file (i.e. otu_table from amptk cluster)
                 -o, --out           Base name for output file. Default: amptk-taxonomy.<method>.txt
                 -d, --db            Select Pre-installed database [ITS1, ITS2, ITS, 16S, LSU, COI]. Default: ITS2
                 -m, --mapping_file  QIIME-like mapping file
                 -t, --taxonomy      Taxonomy calculated elsewhere. 2 Column file.
                 --method            Taxonomy method. Default: hybrid [utax, sintax, usearch, hybrid, rdp, blast]
                 --add2db            Add FASTA files to DB on the fly.
                 --fasta_db          Alternative database of fasta sequenes to use for global alignment.
                 --utax_db           UTAX formatted database. Default: ITS2.udb [See configured DB's below]
                 --utax_cutoff       UTAX confidence value threshold. Default: 0.8 [0 to 0.9]
                 --usearch_db        USEARCH formatted database. Default: USEARCH.udb
                 --usearch_cutoff    USEARCH threshold percent identity. Default 0.7
                 --sintax_cutoff     SINTAX confidence value threshold. Default: 0.8 [0 to 0.9]
                 -r, --rdp           Path to RDP Classifier. Required if --method rdp
                 --rdp_db            RDP Classifer DB set. [fungalits_unite, fungalits_warcup. fungallsu, 16srrna]  
                 --rdp_cutoff        RDP Classifer confidence value threshold. Default: 0.8 [0 to 1.0]
                 --local_blast       Local Blast database (full path) Default: NCBI remote nt database   
                 --tax_filter        Remove OTUs from OTU table that do not match filter, i.e. Fungi to keep only fungi.
                 -u, --usearch       USEARCH executable. Default: usearch9
                 --debug             Keep intermediate files

    Databases Configured: 
    DB_name         DB_type   FASTA originated from                               Fwd Primer   Rev Primer   Records  
    16S.udb         utax      rdp_16s_v16s_sp.fa                                  None         None         13143    
    COI.udb         usearch   arthropods_chordates.bold.fa                        None         None         423813   
    COI_UTAX.udb    utax      coi.genus4utax.fa                                   None         None         60000    
    ITS.udb         usearch   UNITE_public_28.06.2017.fasta                       ITS1-F       ITS4         592807   
    ITS1_UTAX.udb   utax      sh_general_release_dynamic_s_28.06.2017_dev.fasta   ITS1-F       ITS2         57343    
    ITS2_UTAX.udb   utax      sh_general_release_dynamic_s_28.06.2017_dev.fasta   fITS7        ITS4         56015    
    ITS_UTAX.udb    utax      sh_general_release_dynamic_28.06.2017_dev.fasta     ITS1-F       ITS4         30437    
    LSU.udb         usearch   current_Fungi_unaligned.fa                          None         None         92378     
    