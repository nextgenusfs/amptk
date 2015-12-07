#!/usr/bin/env python

import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if len(sys.argv) < 2:
    print "Usage: %s input.fastq BC_name output.fastq\nExample: %s test.demux.fq BC_27 test.output.fq" % (sys.argv[0], sys.argv[0])
    sys.exit(1)


def subtract_sample(file, output):
    global remove_count
    subtract = ";barcodelabel="  + sys.argv[2] + ";"
    with open(output, 'w') as output:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            if title.endswith(subtract):
                remove_count += 1
                continue
            else:
                output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

remove_count = 0
subtract_sample(sys.argv[1], sys.argv[3])
      
print("Removed %i reads that contained %s in the header" % (remove_count, sys.argv[2]))


