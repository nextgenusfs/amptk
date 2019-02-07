#!/usr/bin/env python

from __future__ import print_function
import sys, argparse, os, inspect, itertools, multiprocessing
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

def countBarcodes(file):
    #now loop through data and find barcoded samples, counting each.....
    BarcodeCount = {}
    with open(file, 'r') as input:
        header = itertools.islice(input, 0, None, 4)
        for line in header:
            ID = line.split("=")[-1].split(";")[0]
            if ID not in BarcodeCount:
                BarcodeCount[ID] = 1
            else:
                BarcodeCount[ID] += 1
    return BarcodeCount

def filter_sample(file, output):
    global keep_count, total_count
    with open(output, 'w') as out:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            total_count += 1
            sample = title.split('=',1)[1].split(';')[0]
            if not sample in keep_list:
                keep_count += 1
                if args.format == 'fastq':
                    out.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                elif args.format == 'fasta':
                    out.write(">%s\n%s\n" % (title, seq))

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-remove_samples.py',
		description='''Script parses AMPtk de-multiplexed FASTQ file and keeps those sequences with barocde names in list ''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)
	parser.add_argument('-i','--input', required=True, help='Input AMPtk demux FASTQ')
	parser.add_argument('-l','--list', nargs='+', help='Input list of (BC) names to remove')
	parser.add_argument('-t','--threshold', type=int, help='Keep samples with more reads than threshold')
	parser.add_argument('-f','--file', help='File containing list of names to remove')
	parser.add_argument('-o','--out', required=True, help='Output name')
	parser.add_argument('--format', default='fastq', choices=['fastq','fasta'], help='format of output file')
	args=parser.parse_args(args)

	#check if input compressed, incompress if it is
	if args.input.endswith('.gz'):
		SeqIn = args.input.replace('.gz', '')
		amptklib.Funzip(args.input, SeqIn, multiprocessing.cpu_count())
	else:
		SeqIn = args.input

	remove = []
	if args.threshold:
		print("Finding samples with less than %i reads" % args.threshold)
		BC_counts = countBarcodes(SeqIn)
		for k,v in list(BC_counts.items()):
			if int(v) <= args.threshold:
				if not k in remove:
					remove.append(k)
		print("Removing samples: %s" % ','.join(remove))
	
	if args.file:   
		#load in list of sample names to keep
		with open(args.file, 'r') as input:
			lines = [line.rstrip('\n') for line in input]
		print("Removing samples from file: %s" % ','.join(lines))
		remove = remove + lines

	if args.list:
		lines = args.list
		print("Removing samples from list: %s" % ','.join(lines))
		remove = remove + lines

	#make sure it is a set, faster lookup
	keep_list = set(remove)
	count = len(keep_list)

	#now run filtering 
	keep_count = 0
	total_count = 0

	#rename to base
	if args.out.endswith('.gz'):
		outfile = args.out.replace('.gz', '')
	else:
		outfile = args.out
	#run filtering
	filter_sample(SeqIn, outfile)
	#compress and clean
	if args.out.endswith('.gz'): #compress in place
		amptklib.Fzip_inplace(outfile)
	if args.input.endswith('.gz'):
		amptklib.removefile(SeqIn)

   
	print("Removed %i samples" % count)
	print("Kept %i reads out of %i total reads" % (keep_count, total_count))

if __name__ == "__main__":
	main(args)
