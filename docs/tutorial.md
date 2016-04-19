###UFITS walkthrough

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
As you can see, FASTQ files contain more information than FASTA files.  So we can use the quality scores in the FASTQ file to quality trim the data so that we are removing sequences that have a higher probability of being incorrect.  In NGS sequencing, bases near the 3' end of the sequence are typically of lower quality and thus many quality trimming algorithms start at the 3' end of the molecule and work backwards.  There are many ways to quality trim FASTQ files: average quality score, sliding window quality scores, etc.  In UFITS, we use [Expected Errors](http://drive5.com/usearch/manual/exp_errs.html) trimming, which takes into account the accumulation of expected error across the entire read and throws out reads where the error is greater than a threshold.  This ensures that high quality reads are used to generate OTU clusters and helps to reduce the number of spurious OTUs.

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

####Pre-processing data using UFITS
So you can now see that processing data from different instruments requires a slightly differen strategy, in the case of Illumina data the samples are already split into two files (forward and reverse reads) for each sample, while for Ion/454 the data are combined into a single file that has a unique barcode at the start of the read.  

