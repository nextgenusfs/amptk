#!/usr/bin/env python
import sys

#script to reindex reads containing a ';' such as ';barcodelabel=;'

def fastqreindex(input):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    count = 1
    with open(input, 'rU') as fastq:
        for title, sequence, qual in FastqGeneralIterator(fastq):
            cols = title.split(';')
            header = 'R_'+str(count)+';'+cols[1]+';'
            count += 1
            sys.stdout.write("@%s\n%s\n+\n%s\n" % (header, sequence, qual))

fastqreindex(sys.argv[1])