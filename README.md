# UFITS
###USEARCH Fungal ITS Clustering:###

UFITS is a series of scripts to process fungal ITS amplicon data using USEARCH8
___
####Processing Ion Torrent Data:####

From the Torrent Server, analyze the data using the `--disable-all-filters` BaseCaller argument.  This will leave the adapters/key/barcode sequence intact.  The data need to be exported as a FASTQ file, or alternatively use a 3rd party tool to convert the BAM output file to FASTQ (i.e. `bedtools bamtofastq -i <BAM> -fq <FASTQ>`).  You can then de-multiplex the data as follows:

`ufits-process_ion.py --barcodes 1,5,24 data.fastq > data.demux.fq`
>This will find Ion barcodes (1, 5, and 24) and relabel header with that information (barcodelabel=BC_5;). By default, it will look for all 96 Ion Xpress barcodes, specifiy the barcodes you used by a comma separated list. Next the script will find and trim both the forward and reverse primers (default is ITS2 region: fITS7 & ITS4), and then finally will trim or pad with N's to a set length (default: 250 bp).  Trimming to the same length is critcally important for USEARCH to cluster correctly, padding with N's after finding the reverse primer keeps short ITS sequences from being discarded.  These options can be customized using: `--fwd_primer`, `--rev_primer`, `--trim_len`, etc.  Type `-h` for all the available options.

####Processing Illumina MiSeq PE Data:####

Paired-end MiSeq data is typically delivered already de-multiplexed into separate read files, i.e. barcode1_R1.fastq & barcode1_R2.fastq.  You can merge the PE reads by running the following:

`ufits-merge_illumina.py --out BC_1 --for barcode1_R1.fastq --rev barcode1_R2.fastq`
>This will run USEARCH8 `-fastq_mergepairs`, then run `-fastq_join` on the sequences that do not overlap, and finally concatenate the results together to give you a single FASTQ file, BC_1.fq.  This script can also handle .gz FASTQ files as input.

The data from MiSeq does not contain barcode sequences, but we still need to remove primers and trim/pad to a set length.  You can do that as follows:

`ufits-process_illumina.py [options] BC_1.fq > BC_1.demux.fq`
>Same as above, all options can be customized: `--fwd_primer`, `rev_primer`, `--trim_len`, etc.  Type `-h` for all the available options.  Default is for the file name to be incorporated into the read header name, like barcodelabel=BC_1.

You will need to run this for each barcode file in your Illumina dataset and then concatenate the demuxed files together.  One way to do this for files in the same folder on UNIX is as follows:

```
for file in *.fq; do name=${file%???}; ufits-process_illumina.py $file > $name.demux.fq; done
cat *.demux.fq > illumina.250.fq
```
####OTU Clustering:####

Now the data from either platform (Ion or Illumina) can be clustered by running the following:

`ufits-OTU_cluster.py --out output_250 --maxee 1.0 -f data.demux.fq`
>This will run USEARCH8 `-fastq_filter`, then `-derep_fulllength`, then `-sortbysize`, and finally `-cluster_otus`.  You can also optionally run UCHIME Reference filtering by adding the `--uchime_ref ITS2` option or change the default clustering radius (97%) by passing the `--pct_otu` option. Another option is to process a spike-in control or mock community, you can specify a barcode name for the mock community by passing in `--mock BC_5` which will run some additional steps and report stats of the run to STDOUT.  Type `-h` for all the available options.

An example of a run with EE=0.5, UCHIME ITS2 filtering, and mock community with BC_5:

`ufits-OTU_cluster.py -o test -e 0.5 --uchime_ref ITS2 --mock BC_5 -f data.demux.fq`

####OTU Table Filtering####

The data may need some additional filtering if you included a spike-in control mock community.  The advantage is that you know what should be in the spike-in control barcode sample, thus you can modify USEARCH8 clustering parameters that give you reasonable results.  If you need to trim your OTU table by some threshold, i.e. several OTUs at low abundance are showing up in your spike-in control sample that represent contamination or sequence error - you can set a threshold and filter the OTU table. This is done with the following script:

`ufits-mock_filter.py --trim_data -i test.otu_table.txt -b BC_27`
>This script will read the OTU table `-i`, count the number of OTUs in the barcode specified by the `-b` parameter and give you some basic stats to STDOUT.  It will then ask for a value to threshold trim the data, if you would type in a value of 2, then 2 will be subtracted from every column and a new OTU table will be saved to file ending in `.filtered.out_table.txt` as well as a new OTU fasta file `filtered.otus.fa`.


####Assign Taxonomy:####

Coming soon....

####Dependencies####
python, biopython, USEARCH8 (accessible in PATH; alternatively you can pass in the variable `-u /path/to/usearch8` to scripts requiring USEARCH8).

