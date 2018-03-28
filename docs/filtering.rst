
.. _filtering:

AMPtk OTU Table Filtering
================
An NGS sequencing artefact where reads are assigned to the wrong barcode sample has been reported several times in the literature.  It's been referred to as "index-hopping", "barcode crossover", "contamination", etc, here I refer to this phenomenon as "index-bleed" -> where a small percentage of reads bleed into other samples.  This phenomenon was first reported on Roche 454 platforms and more recently has been reported on Illumina. The mechanism of index-bleed has yet to be determined, however, it seems to happen during the amplification process of NGS sequencing (emulsion PCR on 454/Ion Torrent or cluster generation on Illumina).  Regardless of the mechanism, the impacts of low level index-bleed for downstream community ecology statistics could be large, especially if presence/absence metrics are used. Coupled with countless examples in the literature that show that read abundances in NGS amplicon experiments do not represent biological abundances, it is important to come up with a solution to deal with index-bleed.

Using spike-in control mock communities is a way to measure the sequencing artefacts as well as bioinformatic steps used in a pipeline.  Spike-in controls allow you to see if the number of OTUs generated from a run/software make sense with what you put in.  It is well-known that PCR amplification will bias your sample abundances and is unpredictable in the sense that all metabarcoding amplicons don't amplify with same efficiency in a complex mixture. AMPtk uses spike-in mock communities to measure the degree of index-bleed in a sequencing run and then conservatively applies that threshold to remove read counts that are within the range of index-bleed from an OTU table. The steps are done on an OTU-basis, meaning that low-abundance OTUs are not indiscrimately dropped solely due to the fact that they didn't PCR amplify or sequence well. 

In AMPtk, this process is done using the ``amptk filter`` command, which takes an OTU table and OTUs in FASTA format (i.e. output from any of the amptk clustering commands). 

.. code-block:: none

    Usage:       amptk filter <arguments>
    version:     1.1.1

    Description: Script filters OTU table generated from the `amptk cluster` command and should 
                 be run on all datasets to combat barcode-switching or index-bleed (as high as 
                 2% in MiSeq datasets, ~ 0.3% in Ion PGM datasets).  This script works best when
                 a spike-in control sequence is used, e.g. Synthetic Mock, although a mock is not required.

    Required:    -i, --otu_table     OTU table
                 -f, --fasta         OTU fasta
         
    Optional:    -o, --out           Base name for output files. Default: use input basename
                 -b, --mock_barcode  Name of barcode of mock community (Recommended)
                 -m, --mc            Mock community FASTA file. Required if -b passed. [synmock,mock1,mock2,mock3,other]
                 -c, --calculate     Calculate index-bleed options. Default: all [in,all]
                 -d, --drop          Sample(s) to drop from OTU table. (list, separate by space)
                 --negatives         Negative sample names. (list, separate by space)
                 --ignore            Ignore sample(s) during index-bleed calc (list, separate by space)
         
    Filtering    -n, --normalize     Normalize reads to number of reads per sample [y,n]. Default: y
                 -p, --index_bleed   Filter index bleed between samples (percent). Default: 0.005
                 -t, --threshold     Number to use for establishing read count threshold. Default: max [max,sum,top5,top10,top25]
                 -s, --subtract      Threshold to subtract from all OTUs (any number or auto). Default: 0
                 --delimiter         Delimiter of OTU tables. Default: tsv  [csv, tsv]
                 --min_reads_otu     Minimum number of reads for valid OTU from whole experiment. Default: 2
                 --min_samples_otu   Minimum number of samples for valid OTU from whole experiment. Default: 1
                 --col_order         Column order (separate by space). Default: sort naturally
                 --keep_mock         Keep Spike-in mock community. Default: False
                 --show_stats        Show OTU stats on STDOUT  
                 --debug             Keep intermediate files.
                 -u, --usearch       USEARCH executable. Default: usearch9 


The steps of ``amptk filter`` are:

    1) Maps OTU sequences to those provided from the mock community (``-m, --mc`` argument)

    2) Parses the OTU table, normalizing the read counts for each sample (optional, but recommended)

    3) Next it calculates the number of reads that bleed into the mock community and the number of reads that bleed from the mock community to the rest of the dataset.  The default setting ``-c all`` is desinged for a synthetic mock, if you have biological mock (i.e. real OTUs that might be in your sample) then you can pass the ``-c in`` option to only look at index-bleed into the mock community sample.
    
    4) Then the index-bleed threshold is calculated for each OTU separately based on ``-t, --threshold`` value and read counts less than the calculated threshold are set to 0.
    
    5) The final output then is the filtered OTU table containing actual read counts (normalization is only used for index-bleed filtering).
    
If you do not have a spike-in mock community in your sample, you can still use ``amptk filter`` by providing an index-bleed percentage (``-p, --index_bleed``) which will over-ride the automated calculation.  A value of ``-p 0.005`` or 0.5% is typically able to remove the effects of index-bleed in most MiSeq Illumina datasets.

