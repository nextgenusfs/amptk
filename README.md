# UFITS
###USEARCH Fungal ITS Clustering:###

UFITS is a series of scripts to process fungal ITS amplicon data using USEARCH8.  It can handle Ion Torrent, MiSeq, and 454 data and is cross-platform compatible (works on Windows, Mac, Linux).
___

<img src="https://github.com/nextgenusfs/ufits/blob/master/docs/ufits.png" width="400">


####Installation:####

* [Mac install instructions](docs/mac_install.md)
* [Ubuntu install instructions](docs/ubuntu_install.md)
* [Windows install instuructions](docs/windows_install.md)


####UFITS Wrapper script####

UFITS comes with a wrapper script for ease of use.  On UNIX, you can call it by simply typing `ufits`, while on windows you may need to type `ufits.py` (not needed if you add the `.py` extension in your PATHEXT, directions [here](http://stackoverflow.com/a/13023969/4386003)).

```
$ ufits
Usage:       ufits <command> <arguments>
version:     0.3.3

Description: UFITS is a package of scripts to process fungal ITS amplicon data.  It uses the UPARSE algorithm for clustering
             and thus USEARCH8 is a dependency.
    
Process:     ion         pre-process Ion Torrent data (find barcodes, remove primers, trim/pad)
             illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, remove primers, trim/pad)
             illumina2   pre-process Illumina data from a single file (assumes Ion/454 read structure: <barcode><f_primer>READ)
             454         pre-process Roche 454 (pyrosequencing) data (find barcodes, remove primers, trim/pad)
             select      select reads from de-multiplexed data
             remove      remove reads from de-multiplexed data
             sample      sub-sample (rarify) de-multiplexed reads per sample
             
Clustering:  cluster     cluster OTUs (using UPARSE algorithm)
             filter      OTU table filtering
             taxonomy    Assign taxonomy to OTUs

Utilities:   summarize   Summarize Taxonomy (create OTU-like tables and/or stacked bar graphs for each level of taxonomy)
             funguild    Run FUNGuild (annotate OTUs with ecological information) 
             meta        pivot OTU table and append to meta data
             heatmap     Create heatmap from OTU table
             SRA         De-multiplex data and create meta data for NCBI SRA submission

Setup:       install     Download/install pre-formatted taxonomy DB (UNITE DB formatted for UFITS). Only need to run once.
             database    Format Reference Databases for Taxonomy
```

And then by calling one of the commands, you get a help menu for each:

```
$ ufits cluster
Usage:       ufits cluster <arguments>
version:     0.3.3

Description: Script is a "wrapper" for the UPARSE algorithm.  Modifications include support for a mock spike-in
             community.  FASTQ quality trimming via expected errors and Dereplication are run in Python which allows
             for the use of datasets larger than 4GB.  Chimera filtering and UNOISE are also options.
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             -o, --out           Output base name. Default: out
             -e, --maxee         Expected error quality trimming. Default: 1.0
             -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
             -m, --minsize       Minimum size to keep (singleton filter). Default: 2
             -l, --length        Length to trim reads. Default 250
             --uchime_ref        Run Chimera filtering. Default: off [ITS1, ITS2, Full]
             --map_filtered      Map quality filtered reads back to OTUs. Default: off
             --skip_quality      Skip quality trimming (e.g. reads are already quality trimmed)
             --unoise            Run De-noising pre-clustering (UNOISE). Default: off
             --size_annotations  Append size annotations to OTU names. Default: off
             -u, --usearch       USEARCH executable. Default: usearch8
             --cleanup           Remove intermediate files.
```

####Processing Ion Torrent Data:####

From the Torrent Server, analyze the data using the `--disable-all-filters` BaseCaller argument.  This will leave the adapters/key/barcode sequence intact.  You can download either the unaligned BAM file from the server or a FASTQ file. You can then de-multiplex the data as follows:

```
ufits ion -i data.bam -o data
ufits ion --barcodes 1,5,24 -i data.fastq -o data
ufits ion --barcode_fasta my_barcodes.fa -i data.fastq -o data
```

This will find Ion barcodes (1, 5, and 24) and relabel header with that information (barcodelabel=BC_5;). By default, it will look for all 96 Ion Xpress barcodes, specifiy the barcodes you used by a comma separated list. You can also pass in a fasta file containing your barcode sequences with properly labeled headers. Next the script will find and trim both the forward and reverse primers (default is ITS2 region: fITS7 & ITS4), and then finally will trim or pad with N's to a set length (default: 250 bp).  Trimming to the same length is critcally important for USEARCH to cluster correctly, padding with N's after finding the reverse primer keeps short ITS sequences from being discarded.  These options can be customized using: `--fwd_primer`, `--rev_primer`, `--trim_len`, etc.

####Processing Roche 454 Data:####
Data from 454 instruments has the same read structure as Ion Torrent: <barcode><primer>Read<primer> and thus can be processed very similarly.  You just need to provide either an SFF file, FASTA + QUAL, or FASTQ files and then you need to specify a multi-fasta file containing the barcodes used in the project.  The data will be processed in the same fashion (see above) as the Ion Torrent Data. For example:

```
#SFF input
ufits 454 -i data.sff --barcode_fasta my454barcodes.fa -o 454project

#FASTA/QUAL input
ufits 454 -i data.fa -q data.qual --barcode_fasta my454barcodes.fa -o 454project
```


####Processing Illumina MiSeq PE Data:####

Paired-end MiSeq data is typically delivered already de-multiplexed into separate read files, that have a defined naming structure from Illumina that looks like this: 

```
<sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz
```

You can processes a folder of Illumina data like this:

```
ufits illumina -i folder_name -o miseqData
```

This will find all files ending with '.fastq.gz' in the input folder, gunzip the files, and then sequentially process the paired read files.  First it will run USEARCH8 `-fastq_mergepairs`, however, since some ITS sequences are too long to overlap you can rescue longer sequences by recovering the the non-merged forward reads by passing the `--rescue_forward` argument.  Alternatively, you can only utilize the forward reads (R1), by passing the `--reads forward` argument.  Next the forward and reverse primers are removed and the reads are trimmed/padded to a set length of clustering. Finally, the resulting FASTQ files for each of the processed samples are concatenated together into a file called `"output".demux.fq` that will be used for the next clustering step.  The script will also output a text file called `"output"-filenames.txt` that contains a tab-delimited output of the sample name as well as [i5] and [i7] index sequences that were used.  The script will produce a folder containing the individual de-multiplexed files named from the `-o, --out` argment.

####OTU Clustering:####

Now the data from either platform (Ion, 454, or Illumina) can be clustered by running the following:

```
ufits cluster -i ufits.demux.fq -o ion --uchime_ref ITS2
```

This quality filter the data based on expected errors, then remove duplicated sequences, sort the output by frequency, and finally `usearch -cluster_otus`.  You can also optionally run UCHIME Reference filtering by adding the `--uchime_ref ITS2` option or change the default clustering radius (97%) by passing the `--pct_otu` option. Type `-h` for all the available options.


####OTU Table Filtering####

The data may need some additional filtering if you included a spike-in control mock community.  The advantage is that you know what should be in the spike-in control barcode sample, thus you can modify USEARCH8 clustering parameters that give you reasonable results.  If you need to trim your OTU table by some threshold, i.e. several OTUs at low abundance are showing up in your spike-in control sample that represent contamination or sequence error - you can set a threshold and filter the OTU table. This is done with the following script:

```
ufits filter -i test.otu_table.txt -f test.final.otus.fa -b mock3 --mc my_mock_seqs.fa
```

This will read the OTU table `-i` and the OTUs `-f` from the `ufits cluster` command.  This script will apply an index-bleed filter to clean-up barcode-switching between samples which happens at a rate of ~ 0.2% in Ion Torrent and as much 0.3% in MiSeq data.  The script first normalizes the OTU table to the number of reads in each sample, then (optionally) using the `-b` sample, it will calculate the amount of index-bleed in the OTU table, finally it will loop through each OTU and change values to 0 that are below the `-index_bleed` filter.  Finally, this script will remove the mock spike in control sample from your dataset - as it should not be included in downstream processing, you can keep mock sequences if desired by passing the `--keep_mock` argument.  The output is a filtered OTU table to be used for downstream processing.

If you do not have a mock community spike in, you can still run the index bleed filter by just running the command without a `-b` argument, such as, which will apply a 0.5% filter on the data:

```
ufits filter -i test.otu_table.txt -f test.final.otus.fa -p 0.005
```


####Assign Taxonomy:####

You can assign taxonomy to your OTUs using UFITS, either using UTAX from USEARCH8.1 or using usearch_global.  The databases require some initial setup before you can use the `ufits taxonomy` command.  The following will get you setup with the UNITE database:

```
ufits install --install_unite
```

This commands will download the newest UNITE curated ITS databases.  It will first download the UNITE curated general release, reformat the UNITE headers to be compatible with UTAX classifier training, trim the data for Full length, ITS1, and ITS2 regions, and then finally will train the UTAX classifier with these data.  The script will then download the UNTIE+INSD database, reformat taxonomy in headers and then create a USEARCH database.  The resulting databases are stored in the `DB` folder of the `ufits` directory and are given the names UTAX.udb and USEARCH.udb respectively.

Issuing the `ufits taxonomy` command will inform you which databases have been properly configured as well as usage instructions:

```
$ ufits taxonomy
Usage:       ufits taxonomy <arguments>
version:     0.3.3

Description: Script maps OTUs to taxonomy information and can append to an OTU table (optional).  By default the script
             uses a hybrid approach, e.g. gets taxonomy information from UTAX as well as global alignment hits from the larger
             UNITE-INSD database, and then parses both results to extract the most taxonomy information that it can at 
             'trustable' levels. UTAX results are used if BLAST-like search pct identity is less than 97 pct.  If pct identity
             is greater than 97 pct, the result with most taxonomy levels is retained.
    
Arguments:   -i, --fasta         Input FASTA file (i.e. OTUs from ufits cluster) (Required)
             -o, --out           Base name for output file. Default: ufits-taxonomy.<method>.txt
             -m, --method        Taxonomy method. Default: hybrid [utax, usearch, hybrid, rdp, blast]
             --utax_db           UTAX formatted database. Default: ITS2.udb [See configured DB's below]
             --utax_cutoff       UTAX confidence value threshold. Default: 0.8 [0 to 0.9]
             --usearch_db        USEARCH formatted database. Default: USEARCH.udb
             --usearch_cutoff    USEARCH threshold percent identity. Default 0.7
             -r, --rdp           Path to RDP Classifier. Required if -m rdp
             --rdp_db            RDP Classifer DB set. [fungalits_unite, fungalits_warcup. fungallsu, 16srrna]  
             --rdp_cutoff        RDP Classifer confidence value threshold. Default: 0.8 [0 to 1.0]
             --local_blast       Local Blast database (full path) Default: NCBI remote nt database   
             --append_taxonomy   OTU table to append taxonomy. Default: none
             --only_fungi        Remove non-fungal OTUs from OTU table.
             -u, --usearch       USEARCH executable. Default: usearch8

Databases Configured: 
DB_name       DB_type   FASTA originated from   Fwd Primer   Rev Primer   Records  
FULL.udb      utax      UNITE.utax.fasta        ITS1-F       ITS4         41414    
ITS1.udb      utax      UNITE.utax.fasta        ITS1-F       ITS2         41306    
ITS2.udb      utax      UNITE.utax.fasta        fITS7        ITS4         42176    
USEARCH.udb   usearch   UNITE.usearch.fasta     None         None         534993 
```

And then you can use the `ufits taxonomy` command to assign taxonomy to your OTUs as well as append them to your OTU table as follows:

```
ufits taxonomy -i data.filtered.otus.fa -o output --append_taxonomy data.final.csv
```

####Summarizing the Taxonomy:####

After taxonomy is appended to your OTU table, you can then generate OTU-like tables for each of your samples at all of the levels of taxonomy (i.e. Kingdom, Phylum, Class, Order, Family, Genus).  At the same time, you can create a QIIME-like stacked bar graph from these data.

```
ufits summarize -i data.taxonomy.otu_table.txt -o data-summary
```

The optional `--graphs` argument will create the stacked bar graphs.  You can save in a variety of formats as well as convert the result to precent of total with the `--percent` argument.


####Dependencies####
* Python 2
* Biopython
* USEARCH8 (to use UTAX you will need at least version 8.1.1756)
* natsort
* pandas
* numpy
* matplotlib
* bedtools (only needed if using Ion Torrent BAM file as input)

Python and USEARCH need to accessible in PATH; alternatively you can pass in the variable `-u /path/to/usearch8` to scripts requiring USEARCH8.  

In order to draw a heatmap or stacked bar graph using `ufits.py heatmap` or `ufits summarize` you will need to have the following python libraries installed: `matplotlib, pandas, numpy`.  They can be installed with pip, i.e. `pip install matplotlib pandas numpy`.