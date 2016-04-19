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
The layout of your amplicons will depend on how you designed your experiment.  The read structures of flowgram sequencers (Roche 454 and Ion Torrent) are typically similar, as sequencing is done in one direction.  In this experimental setup, you use fusion primers containing your normal priming site fused to instrument specific adapters, I will use Ion Torrent data as an example.  In order for your amplicon to be sequenced on the Ion platform, the "forward" read must contain the `A-Adapter`, a key signal (TCAG), and optionally a barocde sequence (10-12 bp)
```
Forward primer (Primer A-key):
    5’-CCATCTCATCCCTGCGTGTCTCCGACTCAG-[Barcode]-template-specific-sequence-3’
    
Reverse primer (Primer P1-key):
    5’-CCTCTCTATGGGCAGTCGGTGAT-template-specific-sequence-3’
```