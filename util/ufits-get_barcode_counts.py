#!/usr/bin/env python

import sys, os, itertools, multiprocessing, glob, shutil
from natsort import natsorted
from Bio import SeqIO

def countfastq(input):
    lines = sum(1 for line in open(input))
    count = int(lines) / 4
    return count

def batch_iterator(iterator, batch_size):
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
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
        trunclen = int(sys.argv[3])
        maxee = float(sys.argv[4])
        rec = rec[:trunclen]
        ee = 0
        for bp, Q in enumerate(rec.letter_annotations["phred_quality"]):
            P = 10**(float(-Q)/10)
            ee += P
        if ee <= maxee:
            rec.name = ""
            rec.description = ""
            yield rec

def worker(file):
    name = file.split(".")[0]
    filter_out = name + '.filter.fq'
    with open(filter_out, 'w') as output:
        with open(file, 'rU') as input:
            SeqRecords = SeqIO.parse(file, 'fastq')
            SeqIO.write(MaxEEFilter(SeqRecords), output, 'fastq')

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
    
def filterSeqs(file, lst):
    with open(file, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fastq')
        for rec in SeqRecords:
            bc = rec.id.split("=")[-1].split(";")[0]
            if bc in lst:
                yield rec

countBarcodes(sys.argv[1])
cpus = multiprocessing.cpu_count()
print "----------------------------------"
#split the input FASTQ file into chunks to process
with open(sys.argv[1], 'rU') as input:
    SeqCount = countfastq(sys.argv[1])
    print('{0:,}'.format(SeqCount) + ' records loaded')
    SeqRecords = SeqIO.parse(sys.argv[1], 'fastq')
    chunks = SeqCount / cpus + 1
    print("Using %i cpus to process data" % cpus)
    #divide into chunks, store in tmp file
    pid = os.getpid()
    folder = 'ufits_tmp_' + str(pid)
    if not os.path.exists(folder):
        os.makedirs(folder)
    for i, batch in enumerate(batch_iterator(SeqRecords, chunks)) :
        filename = "chunk_%i.fq" % (i+1)
        tmpout = os.path.join(folder, filename)
        handle = open(tmpout, "w")
        count = SeqIO.write(batch, handle, "fastq")
        handle.close()

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
catDemux = 'concat.filter.fq'
with open(catDemux, 'w') as outfile:
    for filename in glob.glob(os.path.join(folder,'*.filter.fq')):
        if filename == catDemux:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)
shutil.rmtree(folder)
print "----------------------------------"
countBarcodes('concat.filter.fq')
print "----------------------------------"
#make list of those values larger than threshold
threshold = raw_input("Set threshold to remove samples:  ")
keep = []
for k,v in BarcodeCount.items():
    if v >= int(threshold):
        keep.append(k)

#now loop through filtered records and save those to sys.argv[2]
with open(sys.argv[2], 'w') as output:
    SeqIO.write(filterSeqs(catDemux, keep), output, 'fastq')

countBarcodes(sys.argv[2])
print "----------------------------------"
print "Script finished, output in %s" % sys.argv[2]
#finally remove temp file
os.remove(catDemux)

