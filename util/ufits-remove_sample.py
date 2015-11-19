#!/usr/bin/env python

import sys
import re
from Bio import SeqIO


FileName = sys.argv[1]
subtract = ";barcodelabel="  + sys.argv[2] + ";"
OutFile = sys.argv[3]

def subtract_sample(File):
    global remove_count
    with open(File, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fastq')
        for rec in SeqRecords:
            if rec.id.endswith(subtract):
                remove_count += 1
                continue
            else:
                yield rec

remove_count = 0
with open(OutFile, 'w') as output:
    SeqIO.write(subtract_sample(FileName), output, 'fastq')
      
print("Removed %i reads that contained %s in the header" % (remove_count, sys.argv[2]))


