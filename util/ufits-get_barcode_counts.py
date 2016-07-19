#!/usr/bin/env python

import sys, os, itertools, multiprocessing, glob, shutil, argparse
from natsort import natsorted
from Bio import SeqIO


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-get_barcode_counts.py',
    description='''Script loops through demuxed fastq file counting occurances of barcodes, can optionally quality trim and recount.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input demuxed FASTQ')
parser.add_argument('--quality_trim', action='store_true', help='Quality trim data')
parser.add_argument('-e','--maxee', default=1.0, type=float, help='MaxEE Q-trim threshold')
parser.add_argument('-l','--trunclen', default=250, type=int, help='Read truncation length')
parser.add_argument('-o','--out', help='Output for quality trimmed data')
args=parser.parse_args()


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
        trunclen = args.trunclen
        maxee = args.maxee
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

if args.quality_trim and not args.output:
    print "Error, to run quality trimming you must provide -o, --output"
    os._exit(1)

#main start here
cpus = multiprocessing.cpu_count()
print "----------------------------------"
countBarcodes(args.input)
print "----------------------------------"
if args.quality_trim:
    #split the input FASTQ file into chunks to process
    with open(args.input, 'rU') as input:
        SeqCount = countfastq(args.input)
        print('{0:,}'.format(SeqCount) + ' records loaded')
        SeqRecords = SeqIO.parse(input, 'fastq')
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
    catDemux = args.output
    with open(catDemux, 'w') as outfile:
        for filename in glob.glob(os.path.join(folder,'*.filter.fq')):
            if filename == catDemux:
                continue
            with open(filename, 'rU') as readfile:
                shutil.copyfileobj(readfile, outfile)
    shutil.rmtree(folder)
    print "----------------------------------"
    countBarcodes(args.output)
    print "----------------------------------"
    print "Script finished, output in %s" % args.output


