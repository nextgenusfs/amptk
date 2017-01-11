#!/usr/bin/env python

import sys, argparse, os, inspect
from Bio.SeqIO.QualityIO import FastqGeneralIterator
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.amptklib as amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='amptk-remove_samples.py',
    description='''Script parses AMPtk de-multiplexed FASTQ file and keeps those sequences with barocde names in list ''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input AMPtk demux FASTQ')
parser.add_argument('-l','--list', nargs='+', help='Input list of (BC) names to remove')
parser.add_argument('-f','--file', help='File containing list of names to remove')
parser.add_argument('-o','--out', required=True, help='Output name')
parser.add_argument('--format', default='fastq', choices=['fastq','fasta'], help='format of output file')
args=parser.parse_args()

def filter_sample(file, output):
    global keep_count, total_count
    with open(output, 'w') as out:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            total_count += 1
            sample = title.split('=',1)[1].split(';')[0]
            if not sample in keep_list:
                keep_count += 1
                out.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

if not args.list:
    if not args.file:
        print "Error, you must specifiy a list of barcodes or a file containing barcodes"
        os._exit(1)
if not args.file:
    if not args.list:
        print "Error, you must specifiy a list of barcodes or a file containing barcodes"
        os._exit(1)

if args.list and args.file:
    print "Error, you must specifiy either list of barcodes or a file containing barcodes, not both"
    os._exit(1)

if args.file:   
    count = amptklib.line_count(args.file)
    #load in list of sample names to keep
    with open(args.file, 'rU') as input:
        lines = [line.rstrip('\n') for line in input]

if args.list:
    count = len(args.list)
    lines = args.list

#make sure it is a set, faster lookup
keep_list = set(lines)

#now run filtering 
keep_count = 0
total_count = 0
filter_sample(args.input, args.out)
print("Removed %i samples" % count)
print("Kept %i reads out of %i total reads" % (keep_count, total_count))


