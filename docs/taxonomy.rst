
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

    #OTUID  taxonomy
    OTU1    GSL|100.0|KP776995;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Coniochaetales,f:Coniochaetaceae,g:Coniochaeta
    OTU2    US|0.9389|HF674734;k:Fungi,p:Ascomycota
    OTU3    GS|100.0|AY465463;k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Chaetothyriales,f:Herpotrichiellaceae,g:Phialophora
    OTU4    GS|98.0|KX222656;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales,f:Lasiosphaeriaceae
    OTU5    US|1.0000|FM200433;k:Fungi
    OTU6    GS|100.0|KF850373;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales,f:Chaetomiaceae,g:Trichocladium,s:Trichocladium opacum
    OTU7    GDL|100.0|KF660573;k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Pleosporales
    OTU8    GS|100.0|KX611531;k:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,f:Helotiaceae,g:Varicosporium,s:Varicosporium elodeae
    OTU9    GSL|100.0|JN655624;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Hypocreales,g:Ilyonectria
    OTU10   SS|0.9000|KT197416;k:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales

Taxonomy information follows the ; in the string. The pre-taxonomy identifiers are coded to mean the following:

.. code-block:: none

    #OTUID  taxonomy
    OTU1    Method|Pident/Score|Ref_DB ID;

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

These databases were created from Unite v8.0, first downloading two databases from the UNITE website.  First the General FASTA release of the DB `here <https://unite.ut.ee/sh_files/sh_general_release_28.06.2017.zip>`_, and `here <https://unite.ut.ee/sh_files/sh_general_release_s_28.06.2017.zip>`_.  Then the Full UNITE+INSD database `here <https://unite.ut.ee/sh_files/UNITE_public_28.06.2017.fasta.zip>`_.  For the general FASTA releases, the 'developer' fasta files are used. The taxonomy information is then reformated and databases produced as follows:

.. code-block:: none

    #Create full length ITS USEARCH Database, convert taxonomy, and create USEARCH database
    amptk database -i UNITE_public_all_02.02.2019.fasta -f ITS1-F -r ITS4 \
        --primer_required none -o ITS --create_db usearch --install --source UNITE:8.0
        
    #create SINTAX database
    amptk database -i sh_general_release_dynamic_all_02.02.2019_dev.fasta \
        -o ITS_SINTAX --create_db utax -f ITS1-F -r ITS4 --derep_fulllength \
        --install --source UNITE:8.0 --primer_required none    

    #Create UTAX Databases
    amptk database -i sh_general_release_dynamic_all_02.02.2019_dev.fasta  \
        -o ITS_UTAX --create_db utax -f ITS1-F -r ITS4 \
        --derep_fulllength --install --source UNITE:8.0 --primer_required none
        
    amptk database -i sh_general_release_dynamic_all_02.02.2019_dev.fasta \
        -o ITS1_UTAX -f ITS1-F -r ITS2 --primer_required rev --derep_fulllength \
        --create_db utax --install --subsample 65000 --source UNITE:8.0
        
    amptk database -i sh_general_release_dynamic_all_02.02.2019_dev.fasta \
        -o ITS2_UTAX --create_db utax -f fITS7 -r ITS4 --derep_fulllength \
         --install --source UNITE:8.0 --primer_required for

**Arthropod/Chordate mtCOI DB**

These data were pulled from the `BOLDv4 database <http://v4.boldsystems.org>`_  Since most studies using mtCOI regions are interested in identification of insects in the diets of animals, the BOLD database was queried as follows.  All Chordata sequences were downloaded by querying the `BIN database using the search term Chordata <http://v4.boldsystems.org/index.php/Public_BINSearch?query=Chordata&searchBIN=Search+BINs>`_.  Similarly, the Arthropods were searched by querying the `BIN databases using the search term Arthropoda <http://v4.boldsystems.org/index.php/Public_BINSearch?query=Arthropoda&searchBIN=Search+BINs>`_.  All data was then downloaded as TSV output.

The TSV output files (~ 6GB) where then each formatted using the following method, which reformats the taxonomy information and pulls sequences that are annotated in BINS and then clusters sequences in each bin to 99%.

Since it can literally take days to download the arthropod dataset, if you'd like to experiment with the data you can get a copy here: `chordates <https://osf.io/9bh2f/download?version=1>`_ and `arthropods <https://osf.io/aqrey/download?version=1>`_.

.. code-block:: none

    #reformat taxonomy
    bold2utax.py -i arthropods.bold.02092019.txt -o chordates --cluster 99 --drop_suppressed
    bold2utax.py -i arthropods.bold.02092019.txt -o arthropods --cluster 99 --drop_suppressed

    #combine datasets for usearch
    cat arthropods.bold-reformated.fa chordates.bold-reformated.fa > arth-chord.bold-reformated.fasta
    
    #generate global alignment database
    amptk database -i arth-chord.bold.reformated.fasta -f LCO1490 -r mlCOIintR --primer_required none \
        --derep_fulllength --format off --primer_mismatch 4 -o COI --min_len 200 --create_db usearch \
        --install --source BOLD:20190219

The second set of output files from `bold2utax.py` are named with `.BIN-consensus.fa` which are the result of 99% clustering for each BIN. We will combine those for the two datasets and then use those data to generate the SINTAX and UTAX databases.

.. code-block:: none

    #combine datasets
    cat arthropods.BIN-consensus.fa chordates.BIN-consensus.fa > arth-chord.bold.BIN-consensus.fasta
 
    #generate SINTAX database
    amptk database -i arth-chord.bold.BIN-consensus.fasta -f LCO1490 -r mlCOIintR --primer_required none \
        --derep_fulllength --format off --primer_mismatch 4 -o COI_SINTAX --min_len 200 --create_db sintax \
        --install --source BOLD:20190219
        
    #generate UTAX database, need to subsample for memory issues with 32 bit usearch and we require rev primer match here
    amptk database -i arth-chord.bold.BIN-consensus.fasta -f LCO1490 -r mlCOIintR --primer_required rev \
        --derep_fulllength --format off --subsample 00000 --primer_mismatch 4 -o COI_UTAX --min_len 200 \
        --create_db utax --install --source BOLD:20190219

**LSU database**

The fungal 28S database (LSU) was downloaded from `RDP <http://rdp.cme.msu.edu/download/current_Fungi_unaligned.fa.gz>`_.  The sequences were then converted into AMPtk databases as follows:

.. code-block:: none

    amptk database -i RDP_v8.0_fungi.fa -o LSU --format rdp2utax --primer_required none \
        --skip_trimming --create_db usearch --derep_fulllength --install --source RDP:8

    amptk database -i RDP_v8.0_fungi.fa -o LSU_SINTAX --format rdp2utax --primer_required none \
        --skip_trimming --create_db sintax --derep_fulllength --install --source RDP:8

    amptk database -i RDP_v8.0_fungi.fa -o LSU_UTAX --format rdp2utax --primer_required none \
        --skip_trimming --create_db utax --derep_fulllength --install --source RDP:8 --subsample 45000
        
            
To generate a training set for UTAX, the sequences were first dereplicated, and clustered at 97% to get representative sequences for training.  This training set was then converted to a UTAX database:

.. code-block:: none

    amptk database -i fungi.trimmed.fa -o LSU_UTAX --format off \
        --skip_trimming --create_db utax --keep_all --install

**16S database**
This is downloaded from `R. Edgar's website <http://drive5.com/utax/data/rdp_v16.tar.gz>`_ and then formatted for AMPtk.  Note there is room for substantial improvement here, I just don't typically work on 16S - so please let me know if you want some suggestions on what to do here.  Here I reformatted the "domain" taxonomy level to "kingdom" for simplicity (even though I know it is taxonomically incorrect).

.. code-block:: none

    amptk database -i rdp_16s_v16_sp.kingdom.fa -o 16S --format off --create_db usearch \
        --skip_trimming --install --primer_required none --derep_fulllength
        
    amptk database -i rdp_16s_v16_sp.kingdom.fa -o 16S_SINTAX --format off --create_db sintax \
        -f 515FB -r 806RB --install --primer_required for --derep_fulllength
        
    amptk database -i rdp_16s_v16_sp.kingdom.fa -o 16S_UTAX --format off --create_db sintax \
        -f 515FB -r 806RB --install --primer_required for --derep_fulllength    


Checking Installed Databases
-------------------------------------
A simple ``amptk info`` command will show you all the arguments as well as display which databases have been installed.

.. code-block:: none

    amptk info

    ------------------------------
    Running AMPtk v 1.3.0
    ------------------------------
    Taxonomy Databases Installed:
    ------------------------------
     DB_name         DB_type              FASTA originated from                Fwd Primer Rev Primer Records     Date   
     ITS.udb         usearch                   UNITE_public_01.12.2017.fasta   ITS1-F      ITS4     532025  2018-05-01
     ITS1_UTAX.udb   utax  sh_general_release_dynamic_s_01.12.2017_dev.fasta   ITS1-F      ITS2      57293  2018-05-01
     ITS2_UTAX.udb   utax  sh_general_release_dynamic_s_01.12.2017_dev.fasta    fITS7      ITS4      55962  2018-05-01
     ITS_UTAX.udb    utax    sh_general_release_dynamic_01.12.2017_dev.fasta   ITS1-F      ITS4      30580  2018-05-01
    ------------------------------
   
    