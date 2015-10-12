#!/usr/bin/env python

import sys, os
from Bio import SeqIO

FileName = sys.argv[1]
MinLen = 50
TrimLen = int(sys.argv[2])

input = open(FileName, "rU")
SeqRecords = SeqIO.parse(input, "fasta")
for rec in SeqRecords:
    L = len(rec.seq)
    if L > TrimLen:
        rec.seq = rec.seq[:TrimLen]
    else:
        rec.seq = rec.seq + (TrimLen - L)*'N'
    
    SeqIO.write(rec, sys.stdout, 'fasta')


