#!/usr/bin/env python

import sys, os, itertools, random
from natsort import natsorted
from Bio import SeqIO

def countBarcodes(file):
    global BarcodeCount
    #now loop through data and find barcoded samples, counting each.....
    BarcodeCount = {}
    with open(file, 'rU') as input:
        header = itertools.islice(input, 0, None, 4)
        for line in header:
            ID = line.split("=")[-1].split(";")[0]
            if ID not in BarcodeCount:
                BarcodeCount[ID] = 1
            else:
                BarcodeCount[ID] += 1

    #now let's count the barcodes found and count the number of times they are found.
    barcode_counts = "%20s:  %s" % ('Sample', 'Count')
    for k,v in natsorted(BarcodeCount.items(), key=lambda (k,v): v, reverse=True):
        barcode_counts += "\n%20s:  %s" % (k, str(BarcodeCount[k]))
    print("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

def IndexSeqs(file):
    global SeqIndex
    SeqIndex = SeqIO.index(file, 'fastq')

def filterSeqs(file, lst):
    with open(file, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fastq')
        for rec in SeqRecords:
            if rec.id in lst:
                yield rec

if len(sys.argv) < 2:
    print "Usage: %s input.fastq output.fastq number\nExample: %s test.demux.fq test.rare.fq 10000" % (sys.argv[0], sys.argv[0])
    sys.exit(1)

IndexSeqs(sys.argv[1])
countBarcodes(sys.argv[1])
print "----------------------------------"
print "Now sub-sampling reads down to a max of %s per sample" % sys.argv[3]
Reads = []
for key, value in BarcodeCount.items():
    sample = []
    for rec in SeqIndex:
        ID = rec.split("=")[-1].split(";")[0]
        if key == ID:
            sample.append(rec)
    Reads.append(sample)

Subsample = []
for line in Reads:
    if len(line) > int(sys.argv[3]):
        line = random.sample(line, int(sys.argv[3]))
    Subsample.append(line)

Subsample = [item for sublist in Subsample for item in sublist]

with open(sys.argv[2], 'w') as output:
    SeqIO.write(filterSeqs(sys.argv[1], Subsample), output, 'fastq')

print "----------------------------------"
countBarcodes(sys.argv[2])
print "----------------------------------"
print "Sub-sampling done: %s" % sys.argv[2]
