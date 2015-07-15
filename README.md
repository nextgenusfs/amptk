# ficus
**Fungal ITS Clustering using USEARCH:**
FICUS is a series of scripts to process fungal ITS amplicon data using USEARCH8

**Processing Ion Torrent Data:**

From the Torrent Server, analyze the data using the `--disable-all-filters` BaseCaller argument.  This will leave the adapters/key/barcode sequence intact.  The data need to be exported as a FASTQ file, or alternatively use a 3rd party tool to convert the BAM output file to FASTQ (i.e. bedtools bamtofastq).  You can then de-multiplex the data as follows:

`ficus-process_ion.py --barcodes 1,5,24 data.fastq > data.demux.fq`
>This will find Ion barcodes (1, 5, and 24) and relabel header with that information.  Next it will find and trim both the forward and reverse primers (default is ITS2 region: fITS7 & ITS4), and then finally will trim or pad with N's to a set length (default: 250 bp).  These options can be customized using: `--fwd_primer`, `--rev_primer`, `--trim_len`, etc.  Type `-h` for all the available options.


**Processing Illumina MiSeq PE Data:**

Paired-end MiSeq data is typically delivered already de-multiplexed into separate read files, i.e. barcode1_R1.fastq & barcode1_R2.fastq.  You can merge the PE reads by running the following:

`ficus-merge_illumina.py --out BC_1 barcode1_R1.fastq barcode1_R2.fastq`
>This will run USEARCH8 `-fastq_mergepairs`, then run `-fastq_join` on the sequences that do not overlap, and finally concatenate the results together to give you a single FASTQ file, BC_1.fq.  This script can also handle .gz FASTQ files as input.

The data from MiSeq does not contain barcode sequences, but we still need to remove primers and trim/pad to a set length.  You can do that as follows:

`ficus-process_illumina.py BC_1.fq > BC_1.demux.fq`

You will need to run this for each barcode file in your Illumina dataset and then concatenate the demuxed files together.  One way to do this for files in the same folder on UNIX is as follows:

`for file in *.demux.fq; do name=${file%???}; ficus-process_illumina.py $file > $name.demux.fq; done`

`cat *.deumux.fq > illumina.250.fq`


**OTU Clustering:**

Now the data from either platform (Ion or Illumina) can be clustered by running the following:

`ficus-OTU_cluster.py --out output_250 --maxee 1.0 data.demux.fq`
>This will run USEARCH8 `-fastq_filter`, then `-derep_fulllength`, then `-sortbysize`, and finally `-cluster_otus`.  You can also optionally run UCHIME Reference filtering by adding the `--uchime_ref ITS2` option or change the default clustering radius (97%) by passing the `--pct_otu` option.

**Assign Taxonomy:**

Coming soon....

