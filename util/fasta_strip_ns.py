#!/usr/bin/env python

import sys
from Bio.SeqIO.FastaIO import FastaIterator

def strip_padding(file):
    for record in FastaIterator(open(file)):
        Seq = record.seq.rstrip('N')   
        sys.stdout.write(">%s\n%s\n" % (record.id, Seq))

strip_padding(sys.argv[1])