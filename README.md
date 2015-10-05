# UFITS
###USEARCH Fungal ITS Clustering:###

UFITS is a series of scripts to process fungal ITS amplicon data using USEARCH8
___

![UFITS schematic](docs/ufits.png)


####Installation:####

* [Mac install instructions](docs/mac_install.md)
* [Ubuntu install instructions](docs/ubuntu_install.md)
* [Windows install instuructions](docs/windows_install.md)


####UFITS Wrapper script####

UFITS comes with a wrapper script for ease of use.  On UNIX, you can call it by simply typing `ufits`, while on windows you need to type `ufits.py` (unless you have put the .py extension in your EXPATH).

```
$ ufits.py
Usage:      ufits <command> <arguments>
version:    0.1.1
    
Command:    ion         pre-process Ion Torrent data (find barcodes, trim/pad)
            illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, trim/pad, concatenate)
            cluster     cluster OTUs (using UPARSE algorithm)
            filter      OTU table filtering
            heatmap     Create heatmap from OTU table
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com

```
And then by calling one of the commands, you get a help menu for each:

```
$ ufits.py cluster
Usage:      ufits cluster <arguments>
version:    0.1.1
    
Arguments:  -i, --fastq         Input FASTQ file (Required)
            -o, --out           Output base name. Default: out
            -e, --maxee         Expected error quality trimming. Default: 1.0
            -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
            -m, --minsize       Minimum size to keep (singleton filter). Default: 2
            -l, --length        Length to trim reads. Default 250
            --mock              Name of spike-in mock community. Default: None
            --mc                Mock community FASTA file. Default: ufits_mock3.fa
            --uchime_ref        Run Chimera filtering. Default: off [ITS1, ITS2, Full]
            --map_unfiltered    Map unfiltered reads back to OTUs. Default: off
            --unoise            Run De-noising pre-clustering (UNOISE). Default: off
            --size_annotations  Append size annotations to OTU names. Default: off
            -u, --usearch       USEARCH executable. Default: usearch8
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com

```

####Processing Ion Torrent Data:####

From the Torrent Server, analyze the data using the `--disable-all-filters` BaseCaller argument.  This will leave the adapters/key/barcode sequence intact.  The data need to be exported as a FASTQ file, or alternatively use a 3rd party tool to convert the BAM output file to FASTQ (i.e. `bedtools bamtofastq -i <BAM> -fq <FASTQ>`).  You can then de-multiplex the data as follows:

```
ufits ion --barcodes 1,5,24 -i data.fastq -o data
```

This will find Ion barcodes (1, 5, and 24) and relabel header with that information (barcodelabel=BC_5;). By default, it will look for all 96 Ion Xpress barcodes, specifiy the barcodes you used by a comma separated list. Next the script will find and trim both the forward and reverse primers (default is ITS2 region: fITS7 & ITS4), and then finally will trim or pad with N's to a set length (default: 250 bp).  Trimming to the same length is critcally important for USEARCH to cluster correctly, padding with N's after finding the reverse primer keeps short ITS sequences from being discarded.  These options can be customized using: `--fwd_primer`, `--rev_primer`, `--trim_len`, etc.

####Processing Illumina MiSeq PE Data:####

Paired-end MiSeq data is typically delivered already de-multiplexed into separate read files, that have a defined naming structure from Illumina that looks like this: 

```
<sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz
```

You can processes a folder of Illumina data like this:

```
ufits illumina -i folder_name
```

This will find all files ending with '.fastq.gz' in the input folder, gunzip the files, and then sequentially process the paired read files.  First it will run USEARCH8 `-fastq_mergepairs`, however, since some ITS sequences are too long to overlap you can rescue longer sequences by recovering the the non-merged forward reads.  Next the forward and reverse primers are removed and the reads are trimmed/padded to a set length of clustering. Finally, the resulting FASTQ files for each of the processed samples are concatenated together into a file called `ufits.demux.fq` that will be used for the next clustering step.  The script will also output a text file called `ufits-filenames.txt` that contains a tab-delimited output of the sample name as well as [i5] and [i7] index sequences that were used.

####OTU Clustering:####

Now the data from either platform (Ion or Illumina) can be clustered by running the following:

```
ufits cluster -i ufits.demux.fq -o ion --mock BC_5
```

This will run `usearch -fastq_filter` to filter the data based on expected errors, then remove duplicated sequences, sort the output by frequency, and finally `usearch -cluster_otus`.  You can also optionally run UCHIME Reference filtering by adding the `--uchime_ref ITS2` option or change the default clustering radius (97%) by passing the `--pct_otu` option. Another option is to process a spike-in control or mock community, you can specify a barcode name for the mock community by passing in `--mock BC_5` which will run some additional steps and report stats of the run to STDOUT.  Type `-h` for all the available options.


####OTU Table Filtering####

The data may need some additional filtering if you included a spike-in control mock community.  The advantage is that you know what should be in the spike-in control barcode sample, thus you can modify USEARCH8 clustering parameters that give you reasonable results.  If you need to trim your OTU table by some threshold, i.e. several OTUs at low abundance are showing up in your spike-in control sample that represent contamination or sequence error - you can set a threshold and filter the OTU table. This is done with the following script:

```
ufits filter -i test.otu_table.txt -b BC_27
```

This  will read the OTU table `-i`, count the number of OTUs in the barcode specified by the `-b` parameter and give you some basic stats to STDOUT.  It will then ask for a value to threshold trim the data, if you would type in a value of 2, then 2 will be subtracted from every column and a new OTU table will be saved to file ending in `.filteredX.out_table.txt` as well as a new OTU fasta file `filtered.otus.fa`.  To combat 'barcode switching' or 'index bleed', an additional filter can be run that removes OTU counts that are less than 0.1% of the total for each OTU.  If you used dual indexing on MiSeq and have a lot of indexes that were re-used, you will need to increase this filter to at least 0.5, by passing the argument `-p 0.5` to the script.


####Assign Taxonomy:####

Coming soon....
UTAX, RDP, BLAST

####Dependencies####
python, biopython, USEARCH8 (accessible in PATH; alternatively you can pass in the variable `-u /path/to/usearch8` to scripts requiring USEARCH8).  In order to draw a heatmap using `ufits.py heatmap` you will need to have the following python libraries installed: `matplotlib, pandas, numpy`.  They can be installed with pip, i.e. `pip install matplotlib pandas numpy`.