
.. _pre-processing:

AMPtk Pre-Processing
================================

One of the most critical steps in analyzing amplicon NGS data is in the initial processing steps of a pipeline, i.e. finding barcodes/indexes, removing primers, removing phiX, trimming/padding reads to a desired length, etc.  In AMPtk I refer to these steps as pre-processing steps and the method is slightly different depending on the sequencing platform that you used as well as the experimental design (primer design, sequencing primers, etc).  By being thorough in pre-processing steps you can ensure that only high quality data moves through your pipeline thereby reducing spurious and/or contaminating OTUs. Quality filtering in AMPtk is done at the clustering step and utilizes expected-errors trimming, thus pre-processing steps are mainly focused on removing sequencing artefacts from the amplicon reads, such as primer and barcode sequences.  Additionally, since expected-errors quality trimming is used, the length of the amplicon sequences contributes significantly in this calculation, so care must be taken for amplicons of variable length, such as the internal transcribed spacer (ITS) region used for characterization of fungal communities.

Processing Ion Torrent and 454 Data
-------------------------------------
Roche 454 and Ion Torrent PGM machines are flowgram based sequencers that have similar experimental setup but differ in their nucleotide detection mechanisms. In terms of processing data from these platforms, the pre-processing steps are similar. Flowgram based sequencers sequence the amplicon in one direction, thus generate single-ended (SE) reads. Multiplexing on these instruments is done through incorporation of a unique barcode sequence in-between the system specific adapter sequences and the amplicon-specific primer.  Therefore, reads have the following structure:

.. code-block:: none

    5' - Barcode:For_Primer:amplicon:Rev_primer:adapter - 3'

Typically these read structure is produced via fusion primers that generally look like the following:

.. code-block:: none

    #Ion Torrent Forward primer (Primer A:Key Sequence (TCAG)):
    5'-CCATCTCATCCCTGCGTGTCTCCGACTCAG-[Barcode]--template-specific-primer-3'
    
    #Ion Torrent Reverse primer (Primer P1):
    5'-CCTCTCTATGGGCAGTCGGTGAT--template-specific-primer-3'

This means that we need to pre-process our reads by identifying the unique Barcode at the 5' end of the read which will tell us which multiplexed sample the read belongs to.  Then we need to remove the forward and reverse primers from the reads, keeping in mind that with very long amplicons, we may not be able to get high quality sequence all the way to the 3' end of the read.  However, in a high quality read, we will be able to identify the barcode and find the forward primer -> if we are unable to do this, then the read should be discarded.

The next step in many amplicon pipelines is to truncate the reads to a set length to be used for downstream clustering.  While this approach works when amplicon pools are identical or very close to the same length (such as bacterial 16S and/or LSU sequences), this becomes problematic for amplicons that are variable in length (such as the fungal ITS region). Metabarcoding schemes for fungi use either the ITS1 or ITS2 region, which are each on average of 250 bp, however sequence lengths can vary for reach region from 125 bp up to more than 600 bp.  Therefore, truncating reads at a set length causes one of two problems to occur: i) reads less than the truncated value are discarded resulting in loss of real biological sequences or ii) truncated sequences lose the informative length characteristic as well as information is lost for taxonomy assignment. To address this, AMPtk removes both forward and reverse primers and then truncates reads ONLY if they are longer than a set length (controlled by the ``-l, --trim_len`` argument).  The default value is 300 bp, which is a value where you can get reliable sequence read quality and yet maintain as much information as possible for clustering/taxonomy.

**Summary: for Ion Torrent and 454 platforms, AMPtk pre-processing runs the following steps:**

    1) identify unique barcode (index) sequence

    2) rename sequence header with sample name from barcode

    3) find forward and reverse primers

    4) remove (trim) barcode and primer sequences

    5) if sequence is longer than ``--trim_len``, truncate sequence

These steps are executed with the following commands:

.. code-block:: none

    #process Ion Torrent Data directly from unaligned BAM file
    amptk ion -i mydata.bam -o mydata -f fITS7 -r ITS4 --barcode_fasta barcodes.fa
    
    #process Ion Torrent Data from fastq.gz and trim to 250 bp, specify primer sequences
    amptk ion -i mydata.fasta.gz -o mydata -l 250 -m my_mapping_file.txt \
        -f GTGARTCATCGAATCTTTG -r TCCTCCGCTTATTGATATGC
    
    #process 454 data from SFF output
    amptk 454 --sff MyRun.sff -o mydata -l 270 --barcode_fasta barcodes.fa
    
    #process 454 data from fasta + qual, barcoded on both 5' and 3' ends
    amptk 454 --fasta MyRun.fasta --qual MyRun.qual -o mydata \
        --barcode_fasta fwdbarcodes.fa --reverse_barcode revbarcodes.fa
    
Processing Illumina Data
------------------------------------- 
In contrast to SE reads generated from 454 and Ion Torrent, Illumina data has the added benefit of being able to generated paired-end (PE) reads.  Initially (ca 2013), Illumina data was limited to short read lengths (< 150 bp), however, starting with the Illumina MiSeq platform read lengths quickly increased and currently Illumina can deliver 2 x 300 bp reads, which offers a combined read length longer than Ion Torrent (400-500 bp).  One way to generate illumina compatible amplicons is through a two-step PCR reaction, where you use your primer of interest fused to a general Illumina adapter, a second barcoding reaction then attaches an barcode sequence and an additional adapter.

.. code-block:: none

    Forward Nested Primer Sequence
    5'- ACACTCTTTCCCTACACGACGCTCTTCCGATCT-template-specific-primer -3'

    Reverse Nested Primer Sequence  
    5'- GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-template-specific-primer -3'

The second barcoding reaction then adds a unique barcode [i5] to the 5' adapter and a barcode [i7] to the 3' adapter sequence.  Thus the final construct schematically looks like this:

.. code-block:: none

    5' - Adapter1-[i5]-Adapter2-For_Primer:amplicon:Rev_Primer-Adapter3-[i7]-Adapter4 - 3'

The machine then does 4 different sequencing reactions: 1) Read 1 (sequence from Adapter2 for 300 bp), 2) Index Read 1 (sequences the i5 barcode), 3) Read 2 (sequence in reverse direction from Adapter3, 4) Index Read 2 (sequences the i7 barcode).  The software then strips the adapter sequences, and the reads will then look like this:

.. code-block:: none

    5' -  For_Primer:amplicon:Rev_Primer - 3'

Illumina software then de-multiplexes the reads based on the index sequences and splits each sample into two files that are named as such: ``<sample name>_<barcode>_L<lane number>_R<read number>_<set number>.fastq.gz``

.. code-block:: none

    #Example naming after bcl2fastq
    Sample1_ATCCTTG_L001_R1_001.fastq.gz
    Sample1_TCTGGTA_L001_R2_001.fastq.gz

Unfortunately, Sequencing centers do not distribute data to customers in the same format.  I've seen data in at least 3 formats (there are probably more - let me know if your data is in a format that AMPtk cannot process), thus there are 3 illumina based commands in AMPtk: ``amptk illumina`` processes a folder of demultiplexed PE reads, ``amptk illumina2`` processes data that contain an additional 5' unique index sequence, and ``amptk illumina3`` which processes data as R1, R2, I1 (where index reads are in separate file).  

AMPtk pre-processing of Illumina runs the same steps as in 454/Ion Torrent with some additional steps needed, namely merging of PE reads and filtering phiX spike-in.  The general workflow for Illumina reads is:

    1) Merge PE reads (use USEARCH or VSEARCH)

    2) rename sequence header with sample name
    
    3) filter reads that are phiX (USEARCH)

    4) find forward and reverse primers (pay attention to ``--require_primer`` argument)

    5) remove (trim) primer sequences

    6) if sequence is longer than ``--trim_len``, truncate sequence

Some examples of how to issue these commands:

.. code-block:: none

    #simple folder of PE MiSeq data
    amptk illumina -i miseq_folder/ -o mydata -f ITS1-F -r ITS2 
    
    #same folder, however, keep only full length sequences
    amptk illumina -i miseq_folder/ -o mydata -f ITS1-F -r ITS2 --full_length
    
    #process folder of data, however, custom sequencing primer used
    amptk illumina -i miseq_folder/ -o mydata -f ITS1-F -r ITS2 --require_primer off
    
    #data is in R1, R2, I1 format
    amptk illumina3 --forward data_R1.fastq.gz --reverse data_R2.fastq.gz --index data_I1.fastq.gz \
        -m mapping_file.txt -o mydata --fwd_primer ITS1-F --rev_primer ITS2

Processing SRA Data
------------------------------------- 
Amplicon data from the NCBI Small Read Archive (SRA) is typically provided in a single FASTQ file per sample, i.e. an experiment with 48 samples will have 48 SRA archives, which could be composed of PE reads or SE reads depending on the platform.  Extracting PE Illumina data for such a project could be done like this:  

For example, `BioProject PRJNA400449 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA400449>`_ is composed of 24 samples.  To download this project using `SRApy <https://github.com/kdmurray91/SRApy>`_ you could run the following:

.. code-block:: none
    
    #download all SRA runs from the BioProject PRJNA400449
    get-project-sras.py -e myemail@address.edu -d output_folder -p 400449 -F {name}.sra
    
    #then convert SRA to FASTQ using sra-tools fastq-dump
    cd output_folder; for file in *.sra; do fastq-dump -F --split-files $file; done; cd ..

Since these files are PE Illumina format, you can use the standard ``amptk illumina`` command:

.. code-block:: none

    amptk illumina -i output_folder -f ITS1-F -r ITS4 --require_primer off -o mydata

On the other hand, Roche 454 data and Ion Torrent data are contained in a single FASTQ file -> notably different than their raw data which is located in a single file with barcode sequences intact.  Typically in SRA downloaded files, the unique barcode sequence is missing, therefore you can get these data into AMPtk using the ``amptk SRA`` command which takes a folder of single FASTQ files as input.  For example:

.. code-block:: none
    
    #download all SRA runs from the BioProject PRJNA305905 (Ion Torrent Run)
    get-project-sras.py -e myemail@address.edu -d output_folder -p 305905 -F {name}.sra
    
    #then convert SRA to FASTQ using sra-tools fastq-dump
    cd output_folder; for file in *.sra; do fastq-dump -F $file; done; cd ..

Now you can demultiplex using the ``amptk SRA`` command:

.. code-block:: none

    amptk SRA -i output_folder -f AGGGTGCTCCTCACAGCCCTGTG -r TGTCCCCGCRGCRMATTTCCTG -o mydata


