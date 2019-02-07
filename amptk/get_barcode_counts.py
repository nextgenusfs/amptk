#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import itertools
import multiprocessing
import glob
import shutil
import argparse
import inspect
from natsort import natsorted
from Bio import SeqIO
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

def batch_iterator(iterator, batch_size):
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = next(iterator)
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def MaxEEFilter(records):
    for rec in records:
        rec = rec[:args.trunclen]
        ee = 0
        for bp, Q in enumerate(rec.letter_annotations["phred_quality"]):
            P = 10**(float(-Q)/10)
            ee += P
        if ee <= args.maxee:
            rec.name = ""
            rec.description = ""
            yield rec

def worker(file):
    name = file.split(".")[0]
    filter_out = name + '.filter.fq'
    with open(filter_out, 'w') as output:
        with open(file, 'r') as input:
            SeqRecords = SeqIO.parse(file, 'fastq')
            SeqIO.write(MaxEEFilter(SeqRecords), output, 'fastq')

def countBarcodes(file):
    global BarcodeCount
    #now loop through data and find barcoded samples, counting each.....
    BarcodeCount = {}
    with open(file, 'r') as input:
        header = itertools.islice(input, 0, None, 4)
        for line in header:
            ID = line.split("=",1)[-1].split(";")[0]
            if ID not in BarcodeCount:
                BarcodeCount[ID] = 1
            else:
                BarcodeCount[ID] += 1
            

    #now let's count the barcodes found and count the number of times they are found.
    barcode_counts = "%10s:  %s" % ('Sample', 'Count')
    for k,v in natsorted(list(BarcodeCount.items()), key=lambda k_v: k_v[1], reverse=True):
        barcode_counts += "\n%10s:  %s" % (k, str(BarcodeCount[k]))
    print("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

def getSeqLength(file):
    seqlength = {}
    with open(file, 'r') as input:
        header = itertools.islice(input, 1, None, 4)
        for line in header:
            length = len(line) - 1
            if not length in seqlength:
                seqlength[length] = 1
            else:
                seqlength[length] += 1
    lengthlist = []
    countlist = []
    for k,v in natsorted(list(seqlength.items())):
        lengthlist.append(k)
        countlist.append(v)
    print("Read count: %i" % sum(countlist))
    if len(lengthlist) < 2:
        print("Read Length: %i bp" % lengthlist[0])
    else:
        print("Read length average: %i" % (sum(lengthlist) / len(lengthlist)))
        print("Read length range: %i - %i bp" % (lengthlist[0], lengthlist[-1]))
        
    
def filterSeqs(file, lst):
    with amptklib.gzopen(file, 'r') as input:
        SeqRecords = SeqIO.parse(input, 'fastq')
        for rec in SeqRecords:
            bc = rec.id.split("=")[-1].split(";")[0]
            if bc in lst:
                yield rec

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-get_barcode_counts.py',
		description='''Script loops through demuxed fastq file counting occurances of barcodes, can optionally quality trim and recount.''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)
	parser.add_argument('-i','--input', required=True, help='Input demuxed FASTQ')
	parser.add_argument('--quality_trim', action='store_true', help='Quality trim data')
	parser.add_argument('-e','--maxee', default=1.0, type=float, help='MaxEE Q-trim threshold')
	parser.add_argument('-l','--trunclen', default=250, type=int, help='Read truncation length')
	parser.add_argument('-o','--out', help='Output for quality trimmed data')
	args=parser.parse_args(args)

	if args.quality_trim and not args.out:
		print("Error, to run quality trimming you must provide -o, --output")
		sys.exit(1)

	#main start here
	cpus = multiprocessing.cpu_count()
	print("----------------------------------")
	tmpinput = 'amptk_show.tmp'
	if args.input.endswith('.gz'):
		amptklib.Funzip(args.input, tmpinput, cpus)
	else:
		tmpinput = args.input
	countBarcodes(tmpinput)
	print("----------------------------------")
	getSeqLength(tmpinput)
	print("----------------------------------")
	if args.quality_trim:
		#split the input FASTQ file into chunks to process
		#split fastq file
		SeqCount = amptklib.countfastq(tmpinput)
		pid = os.getpid()
		folder = 'amptk_tmp_' + str(pid)
		amptklib.split_fastq(tmpinput, SeqCount, folder, cpus*2)    
		#now get file list from tmp folder
		file_list = []
		for file in os.listdir(folder):
			if file.endswith(".fq"):
				file = os.path.join(folder, file)
				file_list.append(file)

		p = multiprocessing.Pool(cpus)
		for f in file_list:
			#worker(f)
			p.apply_async(worker, [f])
		p.close()
		p.join()

		#get filtered results
		catDemux = args.out
		with open(catDemux, 'w') as outfile:
			for filename in glob.glob(os.path.join(folder,'*.filter.fq')):
				if filename == catDemux:
					continue
				with open(filename, 'r') as readfile:
					shutil.copyfileobj(readfile, outfile)
		if catDemux.endswith('.gz'):
			amptklib.Fzip_inplace(catDemux)
		shutil.rmtree(folder)
		print("----------------------------------")
		countBarcodes(args.out)
		print("----------------------------------")
		print("Script finished, output in %s" % args.out)

	if args.input.endswith('.gz'):
		amptklib.removefile(tmpinput)

if __name__ == "__main__":
	main(args)
