#!/usr/bin/env python

import sys, argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-remove_samples.py',
    description='''Script parses UFITS de-multiplexed FASTQ file and keeps those sequences with barocde names in list ''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input UFITS demux FASTQ')
parser.add_argument('-l','--list', required=True, help='Input list of (BC) names to keep')
parser.add_argument('-o','--out', required=True, help='Output name')
args=parser.parse_args()

def filter_sample(file, output):
    global keep_count, total_count
    with open(output, 'w') as output:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            total_count += 1
            sample = title.split('barcodelabel=')[1]
            sample = sample[:-1]
            if not sample in keep_list:
                keep_count += 1
                output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

#load in list of sample names to keep
with open(args.list, 'rU') as input:
    lines = [line.rstrip('\n') for line in input]

#make sure it is a set, faster lookup
keep_list = set(lines)

#now run filtering 
keep_count = 0
total_count = 0
filter_sample(args.input, args.out)
      
print("Kept %i reads out of %i total reads" % (keep_count, total_count))


