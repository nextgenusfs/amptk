###AMPtk walkthrough

####Note this just a simulated walkthrough and is not based on the data in the test_data folder.

####The basics
#####What is a FASTA file?
FASTA format is a widely used format to represent a DNA or protein sequence, the format is very simple in that the header or sequence name is on its own line starting with the `>` character.  The sequence is in the following line and each new sequence starts with another `>` and a unique header name.
```
>Sample1 #this is the fasta header
ATGCGGCGCCGCCGCCCGAACGGTTG
>Sample2 #this is second sequence
TTTTGGGGCCCGGGAGGAAGAGAGGG
```

#####What is a FASTQ file?
Unfortunately, FASTQ format is a non-standardized format ([read more here](http://drive5.com/usearch/manual/fastq_files.html) - but does conform to some simple rules.  Each sequence is contained in 4 lines, where the header line starts with `@` character, the 2nd line is the DNA sequence, 3rd line starts with `+` character, and the 4th line contains the per base quality scores.  Quality scores are derived from the sequencing machine and are ASCII characters that represent the estimated probability of an error.  [read more here](http://drive5.com/usearch/manual/quality_score.html).
```
@ZEXM5:00007:00045
TTAAAGGGCCCCATGACATGAATCATCGAT
+
-)-6/66)---%-6;.-----)--------
```
As you can see, FASTQ files contain more information than FASTA files.  So we can use the quality scores in the FASTQ file to quality trim the data so that we are removing sequences that have a higher probability of being incorrect.  In NGS sequencing, bases near the 3' end of the sequence are typically of lower quality and thus many quality trimming algorithms start at the 3' end of the molecule and work backwards.  There are many ways to quality trim FASTQ files: average quality score, sliding window quality scores, etc.  In AMPtk, we use [Expected Errors](http://drive5.com/usearch/manual/exp_errs.html) trimming, which takes into account the accumulation of expected error across the entire read and throws out reads where the error is greater than a threshold.  This ensures that high quality reads are used to generate OTU clusters and helps to reduce the number of spurious OTUs.

#####Amplicon Read layout
The layout of your amplicons will depend on how you designed your experiment.  It is important to understand your read layout/structure so you can properly trim/process the data to get useful information.  Barcoded amplicons need to be de-multiplexed and primers need to be removed.

######Ion Torrent/454
The read structures of flowgram sequencers (Roche 454 and Ion Torrent) are typically similar, as sequencing is done in one direction.  In this experimental setup, you use fusion primers containing your normal priming site fused to instrument specific adapters, I will use Ion Torrent data as an example.  In order for your amplicon to be sequenced on the Ion platform, the "forward" read must contain the `A-Adapter`, a key signal (TCAG), and optionally a barocde sequence (10-12 bp).  The reverse primer contains a `P1-Adapter`, which is responsible for binding to the ISP (Ion Sphere Particle) during the templating reaction.  
```
Forward primer (Primer A-key):
    5’-CCATCTCATCCCTGCGTGTCTCCGACTCAG-[Barcode]--template-specific-primer-3’
    
Reverse primer (Primer P1-key):
    5’-CCTCTCTATGGGCAGTCGGTGAT--template-specific-primer-3’
```
During sequencing with Ion data, the sequencing primer binds to the `A-Adapter` and the sequence then starts with the key signal (TCAG), followed by a barcode tag, your template specific primer, and then finally the amplicon of interest.  Currently Ion can sequence reads ~ 400-500 bp in length, so for some targets you will sequence all the way to the P1 adapter and others you may not get sequence all the way to the end of the read.  But your data that comes off of the machine, will schematically look something like this:
```
5' - Barcode:For_Primer:amplicon:Rev_primer:adapter - 3'
```
######Illumina MiSeq paired-end
The basic structure of MiSeq reads are similar to Ion/454 in the sense that the amplicon needs to have adapters at each end of the molecule.  The key difference (advantage) is that Illumina can sequence the amplicon from both directions and thus generates two reads for each amplicon, a forward and a reverse read.  Currently, the MiSeq platform can sequence 2 x 300 bp reads allowing good coverage and length of most amplicons.  One way to generate illumina compatible amplicons is through a two-step PCR reaction, where you use your primer of interest fused to a general Illumina adapter, a second barcoding reaction then attaches an barcode sequence and an additional adapter.
```
Forward Nested Primer Sequence
5’ACACTCTTTCCCTACACGACGCTCTTCCGATCT-template-specific-primer 3’

Reverse Nested Primer Sequence  
5’ GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-template-specific-primer 3’
```
The second barcoding reaction then adds a unique barcode [i5] to the 5' adapter and a barcode [i7] to the 3' adapter sequence.  Thus the final construct schematically looks like this:
```
5' - Adapter1-[i5]-Adapter2-For_Primer:amplicon:Rev_Primer-Adapter3-[i7]-Adapter4 - 3'
```
The machine then does 4 different sequencing reactions: 1) Read 1 (sequence from Adapter2 for 300 bp), 2) Index Read 1 (sequences the i5 barcode), 3) Read 2 (sequence in reverse direction from Adapter3, 4) Index Read 2 (sequences the i7 barcode).  The software then strips the adapter sequences, and the reads will then look like this:
```
5' -  For_Primer:amplicon:Rev_Primer - 3'
```
Illumina software then de-multiplexes the reads based on the index sequences and splits each sample into two files that are named as such, `<sample name>_<barcode>_L<lane number>_R<read number>_<set number>.fastq.gz`
```
Sample1_ATCCTTG_L001_R1_001.fastq.gz
Sample1_TCTGGTA_L001_R2_001.fastq.gz
``` 

####Pre-processing data using AMPtk
So you can now see that processing data from different instruments requires a slightly different strategy.  But generally, what we want to accomplish in the pre-processing step is to label the read as to which sample it was derived from, remove forward and reverse primers, and then trim/pad the read to a set length.  First we will go through how AMPtk processes Ion Torrent data.

######Ion Torrent Pre-processing
De-multiplexing, primer removal, and trimming/padding is done using the `amptk ion` command.  By typing `amptk ion` into the terminal you will see a help menu with options for the script.
```
amptk ion -i ion.fastq --barcode_fasta barcodes.fa -f GTGARTCATCGAATCTTTG -r TCCTCCGCTTATTGATATGC -o ion
```
Ok, so what is happening?  You see that we give AMPtk two inputs, 1) a FASTQ file from the machine and 2) a barcode fasta file containing the barcodes that we put on our primers. We also tell this script which primers we used. 
So the FASTQ files look something like this:
```
@ZEXM5:00008:00047
TTCTCATTGAACAGTGATCATCGAATCTTTGAACGCACATTTGCGCCCCCCGCTATTACCTCTAGCGGGCATCCTGTTCGAGCGCATTTCAACCCCCCTCAAAGCCCCCCAGCTTGGTGTTGGGGCCCCTACGGCTCTGCCCGTTAAGGCCCCCCCTGAAAAACGAAAGGTGGCGGGGCGCTTCGACTACGGCGTCCCGATCGCCAGTTCAAGGGGCCACTTATACTTCGCCTAGGGGAGGCCTTCCCGGCGATCCAACCCCCCCCGCCTTAAAACAACCACATCTTTAACCCCAAAGGGTTGACCCTCGGAAGTCAGGTAGGAATACCCCGCTGAAACTTAAAGCATATCAACTAAGCGGAGGAATCCACCGACTGGCCCCATAAGAGAGGCCTGAGACTGCCAAGGCACACAGGGGATG
+
;65---)-)-)---/<;66;;;;;7;;;;1;;7;;;;;<CC;?1H@>>><(-----)--)---66666,66660---)------6---)--)-----%---6)--7776/----)-)---3555%---%----)-------)--)-)-)-00/.-%666666,6566,505509<CC6---.)--------)-----)------3---)--)---%-)---)-----)--5)-4222*22.2.4.44*-)-----)-)-------%--)-)---%-5)-)-------)-)-55%--)--)-)----)---)-)-----)--60-)-----%------)--)--)444444--)---)---)-5)-)--)--)----4)---%222.222--)-3----2-225@73434--------%--2
```
The Barcode FASTA file looks like this:
```
>BC_9
TGAGCGGAAC
>BC_10
CTGACCGAAC
>BC_22
TTCGAGACGC
>BC_27
AACCATCCGC
>BC_28
ATCCGGAATC
>BC_33
TTCTCATTGAAC
```
Each sequence in the FASTQ file is processed and we first look for a valid barcode sequence at the start of each read.  So in this example, we have 6 samples or barcode sequences, the script will search each of these barcodes against the start of the sequence in the FASTQ file and look for a match (the sequences must match exactly).  So you can see that `BC_33` is a match against this first sequence.  So the script will first rename the FASTQ header to reflect that it found a valid barcode and then it will remove the barocde sequence, the resulting sequence looks like this:
```
@R_1;barcodelabel=BC_33;
AGTGATCATCGAATCTTTGAACGCACATTTGCGCCCCCCGCTATTACCTCTAGCGGGCATCCTGTTCGAGCGCATTTCAACCCCCCTCAAAGCCCCCCAGCTTGGTGTTGGGGCCCCTACGGCTCTGCCCGTTAAGGCCCCCCCTGAAAAACGAAAGGTGGCGGGGCGCTTCGACTACGGCGTCCCGATCGCCAGTTCAAGGGGCCACTTATACTTCGCCTAGGGGAGGCCTTCCCGGCGATCCAACCCCCCCCGCCTTAAAACAACCACATCTTTAACCCCAAAGGGTTGACCCTCGGAAGTCAGGTAGGAATACCCCGCTGAAACTTAAAGCATATCAACTAAGCGGAGGAATCCACCGACTGGCCCCATAAGAGAGGCCTGAGACTGCCAAGGCACACAGGGGATG
+
--/<;66;;;;;7;;;;1;;7;;;;;<CC;?1H@>>><(-----)--)---66666,66660---)------6---)--)-----%---6)--7776/----)-)---3555%---%----)-------)--)-)-)-00/.-%666666,6566,505509<CC6---.)--------)-----)------3---)--)---%-)---)-----)--5)-4222*22.2.4.44*-)-----)-)-------%--)-)---%-5)-)-------)-)-55%--)--)-)----)---)-)-----)--60-)-----%------)--)--)444444--)---)---)-5)-)--)--)----4)---%222.222--)-3----2-225@73434--------%--2
```
The next step of the script is to find the forward primer sequence and remove the sequence.  Here the script allows 2 mismatches, since our primer design contained a single `A` linker and then we used fITS7 as our forward primer `GTGARTCATCGAATCTTTG`, it is easy to find this sequence and trim it.  The resulting sequence would look like this:
```
@R_1;barcodelabel=BC_33;
AACGCACATTTGCGCCCCCCGCTATTACCTCTAGCGGGCATCCTGTTCGAGCGCATTTCAACCCCCCTCAAAGCCCCCCAGCTTGGTGTTGGGGCCCCTACGGCTCTGCCCGTTAAGGCCCCCCCTGAAAAACGAAAGGTGGCGGGGCGCTTCGACTACGGCGTCCCGATCGCCAGTTCAAGGGGCCACTTATACTTCGCCTAGGGGAGGCCTTCCCGGCGATCCAACCCCCCCCGCCTTAAAACAACCACATCTTTAACCCCAAAGGGTTGACCCTCGGAAGTCAGGTAGGAATACCCCGCTGAAACTTAAA<b>GCATATCAACTAAGCGGAGGA</b>ATCCACCGACTGGCCCCATAAGAGAGGCCTGAGACTGCCAAGGCACACAGGGGATG
+
;7;;;;;<CC;?1H@>>><(-----)--)---66666,66660---)------6---)--)-----%---6)--7776/----)-)---3555%---%----)-------)--)-)-)-00/.-%666666,6566,505509<CC6---.)--------)-----)------3---)--)---%-)---)-----)--5)-4222*22.2.4.44*-)-----)-)-------%--)-)---%-5)-)-------)-)-55%--)--)-)----)---)-)-----)--60-)-----%------)--)--)444444--)---)---)-5)-)--)--)----4)---%222.222--)-3----2-225@73434--------%--2
```

Now the script will look for the reverse primer sequence.  Here we used the ITS4 reverse primer `TCCTCCGCTTATTGATATGC`, so we are looking for the reverse complement of that primer `GCATATCAATAAGCGGAGGA`.  The script will then trim the primer as well as any trailing bases.  The resulting sequence looks like this:
```
@R_1;barcodelabel=BC_33;
AACGCACATTTGCGCCCCCCGCTATTACCTCTAGCGGGCATCCTGTTCGAGCGCATTTCAACCCCCCTCAAAGCCCCCCAGCTTGGTGTTGGGGCCCCTACGGCTCTGCCCGTTAAGGCCCCCCCTGAAAAACGAAAGGTGGCGGGGCGCTTCGACTACGGCGTCCCGATCGCCAGTTCAAGGGGCCACTTATACTTCGCCTAGGGGAGGCCTTCCCGGCGATCCAACCCCCCCCGCCTTAAAACAACCACATCTTTAACCCCAAAGGGTTGACCCTCGGAAGTCAGGTAGGAATACCCCGCTGAAACTTAAA
+
;7;;;;;<CC;?1H@>>><(-----)--)---66666,66660---)------6---)--)-----%---6)--7776/----)-)---3555%---%----)-------)--)-)-)-00/.-%666666,6566,505509<CC6---.)--------)-----)------3---)--)---%-)---)-----)--5)-4222*22.2.4.44*-)-----)-)-------%--)-)---%-5)-)-------)-)-55%--)--)-)----)---)-)-----)--60-)-----%------)--)--)
```
The last part of the pre-processing script is to trim/pad the sequence to a set length, this ensures that all sequences are the exact same length which is a requirement for UPARSE clustering as terminal mismatches count (unlike BLAST for example).  Alternatively, one could only use full length amplicon sequences where both forward and reverse primers are found.  The default setting for AMPtk is to trim/pad to 250 bp, sequences that are now shorter than 250 bp get padded with `N's` while those that are longer are trimmed to 250 bp.  You can set this to whatever length you like, however, there are trade offs with quality filtering, length of sequence, and being able to taxonomically classify an OTU.  So an example of what the final output of `amptk ion` will look like is here:
```
@R_1;barcodelabel=BC_33;
AACGCACATTTGCGCCCCCCGCTATTACCTCTAGCGGGCATCCTGTTCGAGCGCATTTCAACCCCCCTCAAAGCCCCCCAGCTTGGTGTTGGGGCCCCTACGGCTCTGCCCGTTAAGGCCCCCCCTGAAAAACGAAAGGTGGCGGGGCGCTTCGACTACGGCGTCCCGATCGCCAGTTCAAGGGGCCACTTATACTTCGCCTAGGGGAGGCCTTCCCGGCGATCCAACCCCCCCCGCCTTAAAACAACCA
+
;7;;;;;<CC;?1H@>>><(-----)--)---66666,66660---)------6---)--)-----%---6)--7776/----)-)---3555%---%----)-------)--)-)-)-00/.-%666666,6566,505509<CC6---.)--------)-----)------3---)--)---%-)---)-----)--5)-4222*22.2.4.44*-)-----)-)-------%--)-)---%-5)-)-
@R_1;barcodelabel=BC_27;
AACGCACATTTGCGCCCCCCGCTATTACCTCTAGCGGGCATCCTGTTCGAGCGCATTTCAACCCCCCTCAAAGCCCCCCAGCTTGGTGTTGGGGCCCCTACGGCTCTGCCCGTTAAGGCCCCCCCTGAAAAACGAAAGGTGGCGGGGCGCTTCGACTACGGCGTCCCGATCGCCGATCCAACCCCCCCATCCCCCGCTGAAACTTAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
;7;;;;;<CC;?1H@>>><(-----)--)---66666,66660---)------6---)--)-----%---6)--7776/----)-)---3555%---%----)-------)--)-)-)-00/.-%666666,6566,505509<CC6---.)--------)-----)------3---)--)---%-)---)-----)--5)-4222*2IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```
######Illumina (MiSeq) Pre-processing
In the case of Illumina data the samples are already split into two files (forward and reverse reads) for each sample, thus we don't need to de-multiplex the data.  However, the paired-end reads need to be merged into a single read - there are a number of ways to do this, AMPtk uses the `-fastq_mergepairs` command of USEARCH8, you can read more [here](http://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html).  The command in AMPtk to process Illumina data looks like this:
```
amptk illumina -i miseq -o illumina -f CTTGGTCATTTAGAGGAAGTAA -r GCTGCGTTCTTTATCGATGC \
    --read_length 250 --rescue_forward
```
So, what is this script doing?  The first step is to merge the forward and reverse reads for each sample - the `-fastq_mergepairs` of USEARCH does this by aligning the reads, determining where they overlap, and merging reads that overlap.  Some reads cannot be merged, this could be for a variety of reasons: 1) the amplicon is longer than what could be sequenced and thus tails don't overlap, 2) the reverse read is of poor quality, or 3) both the reads are of poor quality.  AMPtk can try to rescue reads that do not overlap by using just the forward reads which tend to be of higher quality, this is done with the `--rescue_forward` option.  After reads are merged/rescued, the script then looks for forward and reverse primers, trims them off, and then trims/pads to a set length in the same fashion as Ion data above.

######Can I use full-length reads?
Yes you can, you can pass the `--full_length` argument to both `amptk ion` and `amptk illumina` and only full length reads will be kept.  You will see in the next section how there are trade-offs with virtually every "option" that you choose.  You can elect to use only full length ITS sequences (those that have a valid forward and reverse primer) and then the scripts will not trim/pad the reads to a set length.  In theory, this should be the best way to process the data, however additional bias can be introduced during quality trimming - which is the next step in data processing.

#####Clustering data into OTUs
From here on, data from any of the sequencing platforms is treated identically in AMPtk.  Clustering in AMPtk is done with the `amptk cluster` command, as always the command without any options brings up a help menu.  This command is a wrapper for the UPARSE algorithm.  The first step of `amptk cluster` is to quality filter the data based on accumulation of expected errors (see link above).  Expected errors are calculated for each read, this is done by adding together the probability of error for each base pair of the read.  Therefore a low expected error value equals a high quality read.  So you can easily see why the length of the reads matters a lot for expected errors, as the longer the read, the more expected error you will see.  Here is an example of expected errors for two FASTQ reads:
```
@read22;ee=16;
CGGAAGGATGCAGTGAATCATGAATCTTGAAGCACATTGCGCCCATTAGTATTCTAGTGGGCATCCCTGTTCAGCGCATCAACCCCTCAAAGCGCTTTTTAGTGCTTGGTGTTTGGGGGCTTTCTGCTGCTTGCGCAGGCCCCTGAAAGATAGTGGGCGAGGTCTCCCGTGACCCCTGAGCGTAGAGTTTACACCCTCGCTTTGTGGCTCGCGGCGGGGTTCTCAGCCGTGAAAACCCCCCCCAATTCAC
+
00+/)-)---------366;;56;6660660------)--666,66;A>@@CACC>@><<+000--)---)---66666--)-66/----)---6----%------)-3A5--)-...(---)--------)------3555/A---)-------6/5---)-----)-----55%44-----------)----3)------)-6D0------)----%-)----553------%-------%-)-)-44
@read23;ee=0.38;
TTCCACTTCGCAGTGAATCATCGAATCTTTGAACGCACATTGCACCTACCAGTATTCTGGTAGGTATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCTCTGCTTGGTGTTGGGGCCCTACGCCTCGCGCGTAGGCCCTTAAGACTAGTGGCGGACCTTCTGTGATCCCGAGCGTAGTAATTTTTACCGCTTTGGAGACCTGGAGGCACCGGCCGTTAAACCCCTATTTCTCAAGGTTGACCTCGGATCA
+
CAC=;;>7;;>CCBDCACCBB@BB66;;B=BDACCC@@@D>@=;;0;@G?CCCCC@CCCADCCACC??@@???AC@@@CC>??CC?CG@CD<?;;7;<?@@DECBGB?<;6;;;/;C1;>???@@<?CC@@@CHAHD?DACAC?>;BCCC@CC@CC=?;?CBCCCCCC?CCCCCC>;<<7<BCC5CD@A?CE?F@CCCCACC@CC@CCCACAC=?>ADD?DEE8C>;;1;AACACAC>??DACCCA????
```
An expected error value of 1.0 is equivalent to "the most probably number of errors in a read is zero".  However, in practice this is not the case as errors are estimated from Q scores which are probabilities and are entirely dependent on the sequencing machine outputting accurate Q scores.  Illumina tends to over-estimate it's Q scores and Ion has difficulties sequencing homopolymer regions which results in lower Q scores in these areas.  Thus there is some intrinsic bias based on the sequence of each amplicon and which platform it was sequenced on.  These biases are unavoidable.  Expected errors quality filtering performs well on all platforms and greatly reduces the number of spurious OTUs produced by UPARSE in comparison to other more commonly used quality filtering algorithms.

So a command for `amptk cluster` looks like this:
```
amptk cluster -i ion.demux.fq -o ion
```
As mentioned above the first step is to quality filter the reads, this is controlled with the `-e, --maxee` option and is set by default to `1.0`.  You likely don't need to change this.  The next step of UPARSE is to dereplicate the sequences, in other words, finding identical sequences and counting the number of times each sequence is found.  Following dereplication, the script sorts the resulting sequences in order of most abundant to least abundant and removes sequences that are singletons (not found at least twice in the data).  This is again another quality filtering method to remove sequences that are likely to contain errors.  The data is then fed into the `-cluster_otus` command of USEARCH, which starts with the most abundant sequence and does global alignment with every other sequence forming 'clusters' at 97% identity (controlled by the `-p, --pct_otu` option). This program also runs de novo chimera detection while clustering.  After sequences have been clustered into OTUs, we run reference based chimera detection which is based on a curated UNITE database of ITS sequences (using the `--uchime_ref` option).  Finally, the input FASTQ data is mapped to the final OTUs and an OTU table is generated.  The `barcodelabel=sample1` in the FASTQ header is used to determine which sample the sequence was derived from.  So the output from `amptk cluster` is a FASTA file containing OTUs and an OTU table containing the frequencies at which reads from each sample were found in the dataset.  Your final OTUs will look something like this:
```
>OTU_1
AACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTGTCATTAAATTCTCAACTCCCTTTGATTTCTT
CAAAGGTGAGCTTGGATGTTGGAGGCTTGCCGGCTGCAAAGTCGGCTCCTCTGAAATGCATTGACGAAGGGAGTGTGCAT
GATACGGCCTTCGGTGTGATAATGATCGCCGTGGCTGGCTTGCTGTAGCACCTTTGTTTAATCTTTTTCCATTGGGTTGG
AAAAAGCTTG
>OTU_2
AACGCACATTGCACCTACCAGTATTCTGGTAGGTATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCTCTGCTTGGTGTT
GGGGCCCTACGCCTCGCGCGTAGGCCCTTAAGACTAGTGGCGGACCTTCTGTGATCCCGAGCGTAGTAATTTTTACCGCT
TTGGAGACCTGGAGGCACCGGCCGTTAAACCCCTATTTCTCAAGGTTGACCTCGGATCAGGTAGGAATACCCGCTGAACT
TAA
```
And an OTU table is tab-delimited and looks like this:
```
OTUId	BC_28	BC_94	BC_93	BC_27	BC_10	BC_22	BC_48	BC_9	BC_63	BC_38
OTU_1	1834	0	1196	2826	1877	405	0	1903	0	189
OTU_11	0	131	0	0	0	0	0	0	0	0
OTU_26	6	0	2	4	55	2	2	37	0	0
OTU_14	0	0	0	0	0	639	0	1	0	0
OTU_38	2	0	0	1	0	10	0	0	0	0
OTU_25	29	0	74	29	577	55	0	544	0	5
OTU_2	28	1	21	64	94	7	1916	90	115	1
OTU_6	0	98	0	0	0	0	0	0	0	0
OTU_45	1	0	1	1	0	0	28	0	0	0
```
One of the many things that we learned from sequencing mock communities, is that the counts in an OTU table are somewhat meaningless and do not represent the abundances of DNA sequences that were initially in your sample.  The bias is largely driven by the initial PCR reaction on a mixed community.  Either way, we feel the the most conservative approach is to treat the data as presence/absence data as opposed to inferring abundance from the counts in the OTU tables.

#####Filtering OTU table
Despite all of the quality filtering employed to get to this point, there are still some errors that persist in the dataset.  One of these errors is index-bleed or barcode switching, a phenomenon where a barcode switches or bleeds over into another sample where it doesn't belong.  This likely occurs during the templating reaction in Ion Torrent data and seems to occur in Illumina data during the index de-multiplexing step - we have seen as much as 0.3% index bleed in Illumina data and ~ 0.2% index bleed in Ion data.  In order to combat this phenomenon, I engineered a synthetic mock community composed of ITS-like sequences that can be used as a spike-in control for every sequencing run.  These synthetic mock sequences can then be used to measure the amount of index-bleed that is present in your dataset and allows you to set meaningful filtering thresholds to filter out counts that are likely to be derived from experimental error (barcode switching). You use the `amptk filter` command as follows:
```
amptk filter -i ion.otu_table.txt -f ion.final.otus.fa -b mock
```
Ok, now what is it doing?  The script first maps the FASTA sequences to identify which OTUs correspond to the synthetic mock (or actually any mock that you added).  It then removes OTUs that are represented by only a single read (singleton OTUs) which are likely to be spurious.  The script next normalizes the OTU counts against the number of reads in each sample and then calculates the index-bleed from the synthetic mock into other samples as well as the index-bleed from the samples into the mock.  It takes the larger of these two values, rounds up, and then filters the data to remove counts that are less than the calculated index-bleed for each OTU.  Finally, the script removes the mock-spike in sample from the dataset and converts the final filtered OTU table to binary to represent presence/absence.  The script then reports several OTU tables that were created during filtering as well as an updated FASTA file containing the filtered OTUs.

######What if I didn't use a spike-in control?
It's okay, you won't be subject to public humiliation (yet).  You can still use the `amptk filter` command and I would recommend that you do as barcode switching likely happened in your dataset, you just don't have the necessary data to prove it.  You can pass in index-bleed percent to filter like so:
```
amptk filter -i illumina.otu_table.txt -f illumina.final.otus.fa -p 0.005
```
This will apply an index-bleed filter of 0.5% across your OTU table, which should be enough to remove barcode switching from your dataset.

#####Assigning Taxonomy to your OTUs
Ok, almost finished.  Now we have a filtered OTU table and we have the corresponding sequences in FASTA format, now wtf are these damn things!  Traditionally BLAST has been the tool of choice for assigning taxonomy to OTUs, while this can work, BLAST is really not the right tool for the job.  BLAST uses `local alignment` where sequences are chopped up into smaller `words` and then alignments of those smaller words are done - this is useful as it is much faster than `global alignment`, however we want to compare the entire OTU sequence not just a portion of it.  Recently, there have been tools developed that can classify sequences according to a trained dataset, essentially yielding a `p-value` of how likely your OTU sequence matches up agains a reference dataset.  This was first implemented in the RDP Classifier and more recently in USEARCH via the UTAX commands.  Assigning taxonomy in AMPtk can be done using any of these techniques and is done with the `amptk taxonomy` command.  However, the default method is a hybrid approach between `global alignment` and the `utax` classifier.  First, here is the command, then I'll explain whats happening:
```
amptk taxonomy -f ion.mock.otus.fa -i ion.final.binary.csv -d ITS2
```
The script takes your OTU sequences and runs UTAX versus a pre-trained UNITE dataset it then also runs `global alignment` using USEARCH against a larger UNITE dataset.  The results are then processed in this fashion for each OTU: 1) if there is not a global alignment result > 97% it will default to the UTAX result, 2) if global alignment is > 97% the script compares the levels of taxonomy between both results and keeps the one with more information.  So in the `hybrid` method, the script does it's best to find the most taxonomy information it can at trustable levels.  The taxonomy information can then be appended to an OTU table using the `--append_taxonomy` option.  Additionally, non-fungal OTUs can be filtered out of the dataset using the `--only_fungi` option.  So now we have an OTU table that looks like this:
```
OTUId	BC_9	BC_10	BC_22	BC_27	BC_28	BC_38	BC_48	BC_63	BC_93	Taxonomy
OTU_1	1	1	1	1	1	1	0	0	1	UTAX;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Boletales,f:Coniophoraceae,g:Coniophora
OTU_2	1	1	1	1	1	1	1	1	1	UTAX;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales
OTU_4	0	0	1	0	0	1	0	0	0	UTAX;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Cantharellales,f:Hydnaceae,g:Sistotrema,s:Sistotrema brinkmannii
OTU_10	1	1	0	0	1	0	1	1	1	EU552160;k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Diaporthales,f:Togniniaceae
OTU_12	0	0	1	0	0	1	0	0	0	AY590788;k:Fungi,p:Ascomycota
OTU_13	0	0	0	0	0	1	1	0	0	KJ827975;k:Fungi,p:Ascomycota,c:Dothideomycetes
OTU_14	1	0	1	0	0	0	0	0	0	UTAX;k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Hysteriales,f:Gloniaceae,g:Ceno
```

#####Great, now what else can I do with AMPtk?
You can use a few other tools that are built into AMPtk to help you analyze your data.  One of them is to assign functional information using the FUNGuilds database.  This is done with an OTU table that you have appended taxonomy information as follows:
```
amptk funguild -i amptk-taxonomy.otu_table.taxonomy.txt -o funguilds.txt
```
This script will use the Guilds method that is described [here](https://cbs.umn.edu/sites/cbs.umn.edu/files/public/downloads/Nguyenetal2015b.pdf).

You can also summarize the taxonomy in your OTUs by using the `amptk summarize` command, like so:
```
amptk summarize -i amptk-taxonomy.otu_table.taxonomy.txt -o summary --graphs --format pdf
```
This script will traverse the levels of taxonomy in your OTU table and count up the number of times each term is found.  It creates `csv` files of these data and then can optionally make stacked bar graphs from these data with the `--graphs` option.

Downstream community ecology can be done with the OTU table, to analyze the data using the `vegan` package in `Rstats`, it can be beneficial to pivot the OTU table and append it to meta data about your samples.  This can be done with the `amptk meta` command and requires a `csv` meta data file that has sample names in the first column followed by additional columns containing other information, such as treatments.  It can be run like so:
```
amptk meta -i amptk-taxonomy.otu_table.taxonomy.txt -m meta_data.csv -o meta_otu_table.csv
```

Submitting data to SRA is now a lot easier using AMPtk.  I've build an accessory script that can help you get your data into the correct format for deposition into the small read archive (SRA) of NCBI.  The script takes your original data, minimally processes it, and writes an SRA submission file.  Before using the script, you should start a BioProject on NCBI and submit your BioSamples corresponding to the samples in your dataset.  The BioSample worksheet which gets emailed to you as confirmation from NCBI can be fed to the `amptk SRA` script and will autopopulate an SRA submission form for you.  Usage:
```
amptk SRA -i ion.fastq -o sra --barcode_fasta barcodes.fa --biosample BioSample.txt \
    --platform ion --fwd_primer fITS7 --rev_primer ITS4
```

