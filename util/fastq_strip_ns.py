#!/usr/bin/env python

import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def strip_padding(file):
    for title, seq, qual in FastqGeneralIterator(open(file)):
        Seq = seq.rstrip('N')
        Qual = qual[:len(Seq)]
        assert len(Seq) == len(Qual)    
        sys.stdout.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))

strip_padding(sys.argv[1])