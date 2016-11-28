# UFITS
###USEARCH Fungal ITS Clustering:###

UFITS is a series of scripts to process fungal ITS amplicon data using USEARCH8, although it can also be used to process any NGS amplicon data and includes databases setup for analysis of fungal ITS, fungal LSU, bacterial 16S, and insect COI amplicons.  It can handle Ion Torrent, MiSeq, and 454 data and is cross-platform compatible (works on Mac, Linux - and could work on Windows).
___

<img src="https://github.com/nextgenusfs/ufits/blob/master/docs/ufits.png" width="400">


####Installation:####

* [Mac install instructions](docs/mac_install.md)
* [Linux install instructions](docs/ubuntu_install.md)
* [Windows install instuructions](docs/windows_install.md) - Note use on Windows is not recommended.


####UFITS Wrapper script####

UFITS comes with a wrapper script for ease of use.  On UNIX, you can call it by simply typing `ufits`, while on windows you may need to type `ufits.py` (not needed if you add the `.py` extension in your PATHEXT, directions [here](http://stackoverflow.com/a/13023969/4386003)).

```
$ ufits
Usage:       ufits <command> <arguments>
version:     0.6.0

Description: UFITS is a package of scripts to process NGS amplicon data.  It uses the UPARSE algorithm for clustering
             and thus USEARCH is a dependency.
    
Process:     ion         pre-process Ion Torrent data (find barcodes, remove primers, trim/pad)
             illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, remove primers, trim/pad)
             illumina2   pre-process Illumina data from a single file (read structure: <barcode><f_primer>READ<r_primer>)
             454         pre-process Roche 454 (pyrosequencing) data (find barcodes, remove primers, trim/pad)
             show        show number or reads per barcode from de-multiplexed data
             select      select reads (samples) from de-multiplexed data
             remove      remove reads (samples) from de-multiplexed data
             sample      sub-sample (rarify) de-multiplexed reads per sample
             
Clustering:  cluster     cluster OTUs (using UPARSE algorithm)
             dada2       run dada2 denoising algorithm, produces "inferred sequences" (requires R, dada2, ShortRead)
             unoise2     run UNOISE2 denoising algorithm
             cluster_ref closed/open reference based clustering (EXPERIMENTAL)

Utilities:   filter      OTU table filtering
             taxonomy    Assign taxonomy to OTUs
             summarize   Summarize Taxonomy (create OTU-like tables and/or stacked bar graphs for each level of taxonomy)
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
version:     0.6.0

Description: Script is a "wrapper" for the UPARSE algorithm. FASTQ quality trimming via expected 
             errors and Dereplication are run in vsearch if installed otherwise defaults to Python 
             which allows for the use of datasets larger than 4GB.  
             Chimera filtering and UNOISE are also options.
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             -o, --out           Output base name. Default: out
             -e, --maxee         Expected error quality trimming. Default: 1.0
             -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
             -m, --minsize       Minimum size to keep (singleton filter). Default: 2
             --uchime_ref        Run Ref Chimera filtering. Default: off [ITS, LSU, COI, 16S, custom path]
             --map_filtered      Map quality filtered reads back to OTUs. Default: off
             --unoise            Run De-noising pre-clustering (UNOISE). Default: off
             --debug             Keep intermediate files.
             -u, --usearch       USEARCH executable. Default: usearch9
```
####Installing Databases:####
UFITS is pre-configured to deal with amplicons from fungal ITS, fungal LSU, bacterial 16S, and insect COI.  After installation of UFITS, you can download and install these databases using the `ufits install` command (to overwrite existing databases, use the `--force` option.  To install all of the databases, you would type:

```
#install all databases
ufits install -i ITS LSU 16S COI

#install only ITS databases
ufits install -i ITS
```

The resulting databases are stored in the `DB` folder of the `ufits` directory.  These data are used for both Chimera reference filtering and for assigning taxonomy.

####What processing script do I use??####
You need to be familiar with your read structure! For example, do you have barcodes at the 5' end of your amplicons?  Is the primer sequence in your reads or has it been removed by the sequencing software?  I prefer that data be minimally processed - this acts largely as a quality control filter.  For example when using Ion Torrent data (see instructions below) it is beneficial to turn off all filtering of the data on the server, by doing this you ensure that the reads are properly quality filtered and barcodes/primers are used as an additional quality filter.  If you have PE MiSeq data, typically it comes in an already demultiplexed format, i.e. a folder of PE reads for each sample (_R1.fastq.gz, _R2.fastq.gz) - these data can be used directly with the `ufits illumina` command.  A somewhat common (although I think not ideal) approach is to use a custom sequencing primer with Illumina data - and then the output from the sequencer is already demulitplexed and primers are already stripped - in this case, you can use the `ufits illumina` command, but make sure to pass the `--require_primer off` option.   However, I have seen several other formats of data as well - even with Illumina data being in a single file, where then the `ufits illumina2` script is helpful.  Bottom line is if you don't know how your data is structured in terms of barcodes, primers, reads, paired-ends, etc - you need to find out before doing anything else....


####Processing Ion Torrent Data:####

From the Torrent Server, analyze the data using the `--disable-all-filters` BaseCaller argument.  This will leave the adapters/key/barcode sequence intact.  You can download either the unaligned BAM file from the server or a FASTQ file. You can then de-multiplex the data as follows:

```
ufits ion -i data.bam -o data -f fITS7 -r ITS4
ufits ion --barcodes 1,5,24 -i data.fastq -o data -f ITS1-F -r ITS2
ufits ion --barcode_fasta my_barcodes.fa -i data.fastq -o data -f ITS1 -r ITS2
ufits ion -i data.bam -m mapping_file.txt -o data
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
ufits illumina -i folder_name -o miseqData -f fITS7 -r ITS4
ufits illumina -i folder_name -o miseqData -m mapping_file.txt
```

This will find all files ending with '.fastq.gz' in the input folder, gunzip the files, and then sequentially process the paired read files.  First it will run USEARCH8 `-fastq_mergepairs`, however, since some ITS sequences are too long to overlap you can rescue longer sequences by recovering the the non-merged forward reads by passing the `--rescue_forward` argument.  Alternatively, you can only utilize the forward reads (R1), by passing the `--reads forward` argument.  Next the forward and reverse primers are removed and the reads are trimmed/padded to a set length of clustering. Finally, the resulting FASTQ files for each of the processed samples are concatenated together into a file called `miseqData.demux.fq` that will be used for the next clustering step.  The script will also output a text file called `miseqData-filenames.txt` that contains a tab-delimited output of the sample name as well as [i5] and [i7] index sequences that were used.  The script will produce a folder containing the individual de-multiplexed files named from the `-o, --out` argment.

####OTU Clustering:####

The next step is to run `ufits cluster`, which expects de-multiplexed FASTQ data as a single file with `;barcodelabel=Sample_name` in the FASTQ header. If your data is in some other format, you can use other UNIX/Perl/Python scripts to add the `barcodelabel=` to each read and then cluster your data using UFITS.  Note that reads should be [globally trimmed](http://www.drive5.com/usearch/manual/global_trimming.html) and the pre-processing steps in UFITS take steps to ensure high quality data makes it into the clustering algorithm with minimal sequence loss. Now the data from either platform (Ion, 454, or Illumina) can be clustered by running the following:

```
ufits cluster -i ufits.demux.fq -o ion_output
```

This script wil quality filter the data based on expected errors, then remove duplicated sequences (dereplication), sort the output by abundance, and finally cluster using `usearch -cluster_otus` command.  You can also optionally run UCHIME Reference filtering by adding the `--uchime_ref ITS` option or change the default clustering radius (97%) by passing the `--pct_otu` option. Type `-h` for all the available options.

####DADA2 "Clustering":####
Recently there is a new "OTU picking" algorithm for amplicon based datasets called DADA2 that has sensitivity down to single base pairs, see publication [here](http://www.nature.com/nmeth/journal/v13/n7/full/nmeth.3869.html), GitHub [here](https://github.com/benjjneb/dada2).  This algorithm uses a statistical method to infer the original sequence that a read was derived from, foregoing the need to cluster at a set threshold (i.e. 97%).  I've implemented a modified DADA2 pipeline here to work with the current UFITS data structure.  A reminder is that reads for DADA2 must have no N's and have to all length trimmed identically, thus variable length amplicons will be truncated down.  Thus this method is perhaps more suited to something like COI or 16S amplicons.  You can run it as follows:

```
ufits dada2 -i ufits.demux.fq -o dada2_output -l 200
```
The script will quality filter your data, trim for use in DADA2, run DADA2 alogrithm, and then parse the results to output an OTU table and a file containing inferred sequences (OTUs) in fasta format.  These files can be used in all downstream UFITS scripts, i.e. `ufits filter` and `ufits taxonomy`.

####OTU Table Filtering####

The data may need some additional filtering if you included a spike-in control mock community.  The advantage is that you know what should be in the spike-in control barcode sample, thus you can modify USEARCH8 clustering parameters that give you reasonable results.  If you need to trim your OTU table by some threshold, i.e. several OTUs at low abundance are showing up in your spike-in control sample that represent contamination or sequence error - you can set a threshold and filter the OTU table. This is done with the following script:

```
ufits filter -i test.otu_table.txt -f test.final.otus.fa -b mock3 --mc my_mock_seqs.fa
```

This will read the OTU table `-i` and the OTUs `-f` from the `ufits cluster` command.  This script will apply an index-bleed filter to clean-up barcode-switching between samples which happens at a rate of ~ 0.2% in Ion Torrent and as much 0.3% in MiSeq data.  The script first normalizes the OTU table to the number of reads in each sample, then (optionally) using the `-b` sample, it will calculate the amount of index-bleed in the OTU table, finally it will loop through each OTU and change values to 0 that are below the `-index_bleed` filter.  Finally, this script will remove the mock spike in control sample from your dataset - as it should not be included in downstream processing, you can keep mock sequences if desired by passing the `--keep_mock` argument.  The output is a filtered OTU table to be used for downstream processing.

If you do not have a mock community spike in, you can still run the index bleed filter (and you probably should as nearly all NGS data has some degree of barcode switching or index-bleed) by just running the command without a `-b` argument, such as, which will apply a 0.5% filter on the data.  Passing the `-p 0.005` argument will over-ride the calculated index-bleed.

```
ufits filter -i test.otu_table.txt -f test.final.otus.fa -p 0.005
```

####Assign Taxonomy:####

You can assign taxonomy to your OTUs using UFITS, either using UTAX from USEARCH8.1 or using usearch_global.  The databases require some initial setup before you can use the `ufits taxonomy` command.  

Issuing the `ufits taxonomy` command will inform you which databases have been properly configured as well as usage instructions:

```
$ ufits taxonomy
Usage:       ufits taxonomy <arguments>
version:     0.6.0

Description: Script maps OTUs to taxonomy information and can append to an OTU table (optional).  By default the script
             uses a hybrid approach, e.g. gets taxonomy information from SINTAX, UTAX, and global alignment hits from the larger
             UNITE-INSD database, and then parses results to extract the most taxonomy information that it can at 
             'trustable' levels. SINTAX/UTAX results are used if BLAST-like search pct identity is less than 97%.  
             If % identity is greater than 97%, the result with most taxonomy levels is retained.
    
Arguments:   -f, --fasta         Input FASTA file (i.e. OTUs from ufits cluster) (Required)
             -i, --otu_table     Input OTU table file (i.e. otu_table from ufits cluster)
             -o, --out           Base name for output file. Default: ufits-taxonomy.<method>.txt
             -d, --db            Select Pre-installed database [ITS1, ITS2, ITS, 16S, LSU, COI]. Default: ITS2
             -m, --mapping_file  QIIME-like mapping file
             --method            Taxonomy method. Default: hybrid [utax, sintax, usearch, hybrid, rdp, blast]
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
```

And then you can use the `ufits taxonomy` command to assign taxonomy to your OTUs as well as append them to your OTU table as follows:

```
#use of hybrid taxonomy approach
ufits taxonomy -f data.filtered.otus.fa -o output -i data.final.csv -m mapping_file.txt -d ITS2

#filter data to only include OTUs identified to Fungi
ufits taxonomy -f data.filtered.otus.fa -o output -i data.final.csv -m mapping_file.txt -d ITS2 --tax_filter Fungi

#use RDP classifier
ufits taxonomy -f data.filtered.otus.fa -o output -i data.final.csv -m mapping_file.txt --method rdp --rdp_db fungalits_unite -rdp /path/to/classifier.jar
```

####Summarizing the Taxonomy:####

After taxonomy is appended to your OTU table, you can then generate OTU-like tables for each of your samples at all of the levels of taxonomy (i.e. Kingdom, Phylum, Class, Order, Family, Genus).  At the same time, you can create a QIIME-like stacked bar graph from these data.

```
ufits summarize -i data.taxonomy.otu_table.txt -o data-summary
```

The optional `--graphs` argument will create the stacked bar graphs.  You can save in a variety of formats as well as convert the result to precent of total with the `--percent` argument.

####Downstream Processing:####
As of `ufits v0.6.0`, the output from the `ufits taxonomy` command will create a biom file that contains taxonomy compatible with QIIME, [PHINCH](www.phinch.org), [PhyloSeq](https://joey711.github.io/phyloseq/index.html), and [MetaCoMET](http://probes.pw.usda.gov/MetaCoMET/index.php). The script will also output a phylogenetic tree from your OTUs which is required for some downstream analysis.  Moreover, you can pass the `-m, --mapping_file` option to `ufits taxonomy` and all columns will be incorporated as sample metadata.  A mapping file is created automatically for you in the pre-processing steps of ufits, such as `ufits ion` and `ufits illumina`.  You can easily add your metadata to this file using something like excel, and then save as a tab delimited text file.

[QIIME and PhyloSeq import instructions](docs/downstream_processing.md)


####Dependencies####
* Python 2
* Biopython
* USEARCH8 (to use UTAX you will need at least version 8.1.1756)
* natsort
* pandas
* numpy
* matplotlib
* psutil
* bedtools (only needed if using Ion Torrent BAM file as input)
* vsearch (version > 1.9.0, this is optional but will increase speed of UFITS and is required for very large datasets) installed via homebrew installation by default.  Newer tools require vsearch.
* biom-format (to create biom OTU table)
* h5py (for biom)
* R (dada2)
* dada2, ShortRead (these will be automatically installed on first usage of `ufits dada2`

Python and USEARCH need to accessible in PATH; alternatively you can pass in the variable `-u /path/to/usearch8` to scripts requiring USEARCH8.  

In order to draw a heatmap or stacked bar graph using `ufits.py heatmap` or `ufits summarize` you will need to have the following python libraries installed: `matplotlib, pandas, numpy`.  They can be installed with pip, i.e. `pip install matplotlib pandas numpy`.
