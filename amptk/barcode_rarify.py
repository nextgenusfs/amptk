#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import itertools
import random
import argparse
import inspect
import multiprocessing
from natsort import natsorted
from Bio import SeqIO
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

def countBarcodes(file):
    global BarcodeCount
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

    #now let's count the barcodes found and count the number of times they are found.
    barcode_counts = "%20s:  %s" % ('Sample', 'Count')
    for k,v in natsorted(list(BarcodeCount.items()), key=lambda k_v: k_v[1], reverse=True):
        barcode_counts += "\n%20s:  %s" % (k, str(BarcodeCount[k]))
    print("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

def IndexSeqs(file):
    global SeqIndex
    SeqIndex = SeqIO.index(file, 'fastq')

def filterSeqs(file, lst, out):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    with open(out, 'w') as output:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            if title in lst:
               output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-barcode_rarify.py',
		description='''Script to sub-sample reads down to the same number for each sample (barcode)''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)
	parser.add_argument('-i','--input', required=True, help='Input FASTQ')
	parser.add_argument('-n','--num_reads', required=True, type=int, help='Number of reads to rarify down to')
	parser.add_argument('-o','--out', required=True, help='Output name')
	args=parser.parse_args(args)

	#check if input compressed, incompress if it is
	if args.input.endswith('.gz'):
		SeqIn = args.input.replace('.gz', '')
		amptklib.Funzip(args.input, SeqIn, multiprocessing.cpu_count())
	else:
		SeqIn = args.input
	if args.out.endswith('.gz'):
		outfile = args.out.replace('.gz', '')
	else:
		outfile = args.out

	IndexSeqs(SeqIn)
	countBarcodes(SeqIn)
	print("----------------------------------")
	print("Now sub-sampling reads down to a max of %s per sample" % args.num_reads)
	Reads = []
	for key, value in list(BarcodeCount.items()):
		sample = []
		for rec in SeqIndex:
			ID = rec.split("=")[-1].split(";")[0]
			if key == ID:
				sample.append(rec)
		Reads.append(sample)
	print("Finished indexing reads, split up by barcodelabel")
	Subsample = []
	for line in Reads:
		if len(line) > int(args.num_reads):
			line = random.sample(line, int(args.num_reads))
		Subsample.append(line)

	Subsample = [item for sublist in Subsample for item in sublist]

	#convert list to set for faster lookup
	Lookup = set(Subsample)

	print("Finished randomly sampling reads, now writing %i sequences to %s" % (len(Lookup), outfile))
	filterSeqs(SeqIn, Lookup, outfile)
	print("----------------------------------")
	countBarcodes(outfile)
	#compress and clean
	if args.out.endswith('.gz'): #compress in place
		amptklib.Fzip_inplace(outfile)
	if args.input.endswith('.gz'):
		amptklib.removefile(SeqIn)
	print("----------------------------------")
	print("Sub-sampling done: %s" % args.out)

if __name__ == "__main__":
	main(args)
