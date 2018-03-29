#!/usr/bin/env python

#script to reformat green-genes fasta file for UTAX compatibility

from __future__ import print_function
from builtins import str
import sys
from Bio import SeqIO

def dereplicate(input, output):
    seqs = {}
    in_file = open(input, 'rU')
    if input.endswith('.fa') or input.endswith('.fasta'):
        filetype = 'fasta'
    elif input.endswith('.fq') or input.endswith('.fastq'):
        filtype = 'fastq'
    else:
        print("Could not detect file type, must be FASTA or FASTQ")
        sys.exit
    for rec in SeqIO.parse(in_file, filetype):
        sequence = str(rec.seq)
        if sequence not in seqs:
            seqs[sequence] = rec.id
    with open(output, 'w') as out:
        for sequence in seqs:
            out.write('>'+seqs[sequence]+'\n'+sequence+'\n')


dereplicate(sys.argv[1], sys.argv[2])


