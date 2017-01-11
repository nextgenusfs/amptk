#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, logging, subprocess, multiprocessing, glob, itertools, re
from Bio import SeqIO
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.fasta as fasta
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.amptklib as amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-process_illumina_raw.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', dest='fastq', required=True, help='Illumina FASTQ R1 reads')
parser.add_argument('--reverse', required=True, help='Illumina FASTQ R2 reads')
parser.add_argument('--index', required=True, help='Illumina FASTQ index reads')
parser.add_argument('--map_file', help='QIIME-like mapping tool')
parser.add_argument('-o','--out', dest="out", default='illumina_out', help='Base name for output')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--barcode_fasta', default='pgm_barcodes.fa', help='FASTA file containing Barcodes (Names & Sequences)')
parser.add_argument('--reverse_barcode', help='FASTA file containing 3 prime Barocdes')
parser.add_argument('-n','--name_prefix', dest="prefix", default='R_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

def FindBarcode(Seq, BarcodeDict):
    for BarcodeLabel in BarcodeDict.keys():
        Barcode = BarcodeDict[BarcodeLabel]
        if Seq.startswith(Barcode):
            return Barcode, BarcodeLabel
    return "", ""

def TrimRead(record, Ftrim, Rtrim, Name, Count):
    #function to trim a seqrecord and rename
    if Rtrim:
        record = record[Ftrim:Rtrim]
    else:
        record = record[Ftrim:]
    #rename header
    if args.multi == 'False':
        record.id = LabelPrefix + str(Count) + ";barcodelabel=" + Name + ";"
    elif args.multi != 'False':
        if args.multi.endswith('_'):
            args.multi = args.multi.replace('_', '')
        record.id = LabelPrefix + str(Count) + ";barcodelabel=" + args.multi + "_" + Name + ";"
    record.name = ''
    record.description = ''
    return record

def ProcessReads(records):
    global OutCount
    for rec in records:
        #convert to string for processing
        Seq = str(rec.seq)
        
        #look for barcodes
        Barcode, BarcodeLabel = FindBarcode(Seq, Barcodes)
        if Barcode == "": #if not found, try to find with mismatches
            if args.barcode_mismatch > 0:
                hit = [None, None, 0, None, None]
                for k,v in Barcodes.items():
                    alignment = amptklib.fuzzymatch(v, Seq, args.barcode_mismatch)
                    if alignment:
                        if alignment[0] > hit[2]:
                            hit = [k, v, alignment[0], alignment[1], alignment[2]]
                if hit[0] != None:
                    BarcodeLength = hit[4] - hit[3] #might be shorter than actual barcode
                    BarcodeLabel = hit[0]
                    Barcode = hit[1]
                else:
                    continue
            else:
                continue
        else: #barcode was found from dictionary
            BarcodeLength = len(Barcode)

        #now look for primer, if not found, move onto next record
        BestPosFor, BestDiffsFor = primer.BestMatch2(Seq, FwdPrimer, MAX_PRIMER_MISMATCHES)
        if BestPosFor > 0 and BestPosFor <= BarcodeLength+2: #if found will be > 0, and should be found after barcode
            ForTrim = BestPosFor + PL
        else:
            continue
        
        #counter for numbering reads
        OutCount += 1

        #look for reverse primer
        BestPosRev, BestDiffsRev = primer.BestMatch2(Seq, RevPrimer, MAX_PRIMER_MISMATCHES)
        if BestPosRev > 0:  #reverse primer was found    
            #location to trim sequences
            RevTrim = BestPosRev
            
            #determine reverse barcode
            if args.reverse_barcode:
                BCcut = BestPosRev + RL
                CutSeq = Seq[BCcut:]
                if not CutSeq in RevBarcodes:
                    if args.barcode_mismatch > 0:
                        hit = [None, None, 0, None, None]
                        for k,v in RevBarcodes.items():
                            alignment = amptklib.fuzzymatch(k, CutSeq, args.barcode_mismatch)
                            if alignment:
                                if alignment[0] > hit[2]:
                                    hit = [v, k, alignment[0], alignment[1], alignment[2]]
                        if hit[0] != None:
                            BCname = hit[0]
                        else:
                            continue
                    else:
                        continue
                else:
                    BCname = RevBarcodes.get(CutSeq)
                #update name
                BarcodeLabel = BarcodeLabel+'_'+BCname
            
            #trim record
            rec = TrimRead(rec, ForTrim, RevTrim, BarcodeLabel, OutCount)
            
            #check length       
            L = len(rec.seq)
            if L < MinLen:
                continue
            if not args.full_length:
                #now check trim length, pad if necessary
                if L < TrimLen:
                    pad = TrimLen - L
                    Seq = str(rec.seq)
                    Seq = Seq + pad*'N'
                    Qual = rec.letter_annotations["phred_quality"]
                    pad = TrimLen - L
                    add = [40] * pad
                    Qual.extend(add)
                    del rec.letter_annotations["phred_quality"]
                    rec.seq = Seq
                    rec.letter_annotations["phred_quality"] = Qual
                    yield rec
                elif L >= TrimLen:   
                    rec = rec[:TrimLen]
                    yield rec
            else:
                yield rec

        else: #if it is full length, we did not find reverse primer, so drop read
            if not args.full_length:
                #trim record
                rec = TrimRead(rec, ForTrim, False, BarcodeLabel, OutCount)
                #check length
                L = len(rec.seq)
                if L < MinLen: #remove if shorter than minimum length
                    continue
                #truncate down to trim length
                if L >= TrimLen:
                    rec = rec[:TrimLen]        
                    yield rec

def worker(input):
    output = input.split(".",-1)[0] + '.demux.fq'
    with open(output, 'w') as o:
        with open(input, 'rU') as i:
            SeqRecords = SeqIO.parse(i, 'fastq')
            SeqIO.write(ProcessReads(SeqRecords), o, 'fastq')  

args.out = re.sub(r'\W+', '', args.out)

log_name = args.out + '.amptk-demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)
FNULL = open(os.devnull, 'w')
amptklib.setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
amptklib.SystemInfo()
#get version of amptk
version = amptklib.get_version()
amptklib.log.info("%s" % version)

#parse a mapping file or a barcode fasta file, primers, etc get setup
#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
barcode_file = args.out + ".barcodes_used.fa"
if os.path.isfile(barcode_file):
    os.remove(barcode_file)

#check if mapping file passed, use this if present, otherwise use command line arguments
if args.mapping_file:
    if not os.path.isfile(args.mapping_file):
        amptklib.error("Mapping file is not valid: %s" % args.mapping_file)
        sys.exit(1)
    mapdata = amptklib.parseMappingFile(args.mapping_file, barcode_file)
    #forward primer in first item in tuple, reverse in second
    FwdPrimer = mapdata[0]
    RevPrimer = mapdata[1]
    genericmapfile = args.mapping_file
else:
    if args.barcode_fasta:
        if not os.path.isfile(args.barcode_fasta):
            amptklib.error("Mapping file or barcode_fasta is required")
            sys.exit(1)
        else:
            shutil.copyfile(args.barcode_fasta, barcode_file)
    
    #parse primers here so doesn't conflict with mapping primers
    #look up primer db otherwise default to entry
    if args.F_primer in amptklib.primer_db:
        FwdPrimer = amptklib.primer_db.get(args.F_primer)
    else:
        FwdPrimer = args.F_primer
    if args.R_primer in amptklib.primer_db:
        RevPrimer = amptklib.primer_db.get(args.R_primer)
    else:
        RevPrimer = args.R_primer

#setup 
if args.map_file:
    mapdict = amptklib.mapping2dict(args.map_file)
else:
    mapdict = False

print mapdict

#now try to get index reads into dictionary
indexReads = amptklib.barcodes2dict(args.index, mapdict)

print indexReads
sys.exit(1)


if args.barcode_fasta == 'pgm_barcodes.fa':
    amptklib.log.error("You did not specify a --barcode_fasta, it is required this type of data")
    os._exit(1)
if args.reverse:
    FNULL = open(os.devnull, 'w')
    #test for usearch
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        amptklib.log.warning("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    amptklib.log.info("USEARCH version: %s" % usearch_test)
    
    #next run USEARCH8 mergepe
    SeqIn = args.out + '.merged.fq'
    mergelog = args.out + '.mergedPE.log'
    notmerged = args.out + '.notmergedfwd.fq'
    amptklib.log.info("Merging PE Illumina reads with USEARCH")
    amptklib.log.debug("%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastq_truncqual 5 -fastq_maxdiffs 8 -minhsp 12" % (usearch, args.fastq, args.reverse, SeqIn))
    with open(mergelog, 'w') as logfile:
        subprocess.call([usearch, '-fastq_mergepairs', args.fastq, '-reverse', args.reverse, '-fastqout', SeqIn, '-fastq_truncqual', '5','-minhsp', '12','-fastq_maxdiffs', '8', '-fastqout_notmerged_fwd', notmerged], stdout = logfile, stderr = logfile)
    #recover forward reads that were not merged, could be longer sequence and still be ok
    with open(SeqIn, 'a') as output:
        with open(notmerged, 'rU') as input:
            for line in input:
                output.write(line)
else:
    amptklib.log.info("Running AMPtk on forward Illumina reads")
    SeqIn = args.fastq 

RevPrimer = revcomp_lib.RevComp(RevPrimer)
amptklib.log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
barcode_file = args.out + ".barcodes_used.fa"
if os.path.isfile(barcode_file):
    os.remove(barcode_file)

if args.barcode_fasta == 'pgm_barcodes.fa':
    #get script path and barcode file name
    pgm_barcodes = os.path.join(parentdir, 'DB', args.barcode_fasta)
    if args.barcodes == "all":
        shutil.copyfile(pgm_barcodes, barcode_file)
    else:
        bc_list = args.barcodes.split(",")
        inputSeqFile = open(pgm_barcodes, "rU")
        SeqRecords = SeqIO.to_dict(SeqIO.parse(inputSeqFile, "fasta"))
        for rec in bc_list:
            name = "BC_" + rec
            seq = SeqRecords[name].seq
            outputSeqFile = open(barcode_file, "a")
            outputSeqFile.write(">%s\n%s\n" % (name, seq))
        outputSeqFile.close()
        inputSeqFile.close()
else:
    shutil.copyfile(args.barcode_fasta, barcode_file)

#setup barcode dictionary
Barcodes = fasta.ReadSeqsDict(barcode_file)

#setup for looking for reverse barcode
if args.reverse_barcode:
    RevBarcodes = {}
    rev_barcode_file = args.out + '.revbarcodes_used.fa'
    if os.path.isfile(rev_barcode_file):
        os.remove(rev_barcode_file)
    if not os.path.isfile(args.reverse_barcode):
        amptklib.log.info("Reverse barcode is not a valid file, exiting")
        sys.exit(1) 
    shutil.copyfile(args.reverse_barcode, rev_barcode_file)
    #parse and put into dictionary
    with open(rev_barcode_file, 'w') as output:
        with open(args.reverse_barcode, 'rU') as input:
            for rec in SeqIO.parse(input, 'fasta'):
                RevSeq = str(rec.seq.reverse_complement())
                if not RevSeq in RevBarcodes:
                    RevBarcodes[RevSeq] = rec.id
                else:
                    amptklib.log.error("Duplicate reverse barcodes detected, exiting")
                    sys.exit(1)
    
#get number of CPUs to use
if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

#get other values
MAX_PRIMER_MISMATCHES = int(args.primer_mismatch)
LabelPrefix = args.prefix
MinLen = int(args.min_len)
TrimLen = int(args.trim_len)
PL = len(FwdPrimer)
RL = len(RevPrimer)
OutCount = 0

#split the input FASTQ file into chunks to process
with open(SeqIn, 'rU') as input:
    SeqCount = amptklib.countfastq(SeqIn)
    amptklib.log.info('{0:,}'.format(SeqCount) + ' records loaded')
    SeqRecords = SeqIO.parse(SeqIn, 'fastq')
    chunks = SeqCount / cpus + 1
    amptklib.log.info("splitting job over %i cpus, this may still take awhile" % cpus)
    #divide into chunks, store in tmp file
    pid = os.getpid()
    folder = 'amptk_tmp_' + str(pid)
    if not os.path.exists(folder):
        os.makedirs(folder)
    for i, batch in enumerate(amptklib.batch_iterator(SeqRecords, chunks)) :
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
    p.apply_async(worker, [f])
p.close()
p.join()

print "-------------------------------------------------------"
#Now concatenate all of the demuxed files together
amptklib.log.info("Concatenating Demuxed Files")

tmpDemux = args.out + '.tmp.demux.fq'
with open(tmpDemux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(folder,'*.demux.fq')):
        if filename == tmpDemux:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)

#clean up tmp folder
shutil.rmtree(folder)

#last thing is to re-number of reads as it is possible they could have same name from multitprocessor split
catDemux = args.out + '.demux.fq'
amptklib.fastqreindex(tmpDemux, catDemux)
os.remove(tmpDemux)
        
amptklib.log.info("Counting FASTQ Records")
total = amptklib.countfastq(catDemux)
amptklib.log.info('{0:,}'.format(total) + ' reads processed')

#now loop through data and find barcoded samples, counting each.....
BarcodeCount = {}
with open(catDemux, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=")[-1].split(";")[0]
        if ID not in BarcodeCount:
            BarcodeCount[ID] = 1
        else:
            BarcodeCount[ID] += 1

#now let's count the barcodes found and count the number of times they are found.
barcode_counts = "%30s:  %s" % ('Sample', 'Count')
for k,v in natsorted(BarcodeCount.items(), key=lambda (k,v): v, reverse=True):
    barcode_counts += "\n%30s:  %s" % (k, str(BarcodeCount[k]))
amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

#get file size
filesize = os.path.getsize(catDemux)
readablesize = amptklib.convertSize(filesize)
amptklib.log.info("Output file:  %s (%s)" % (catDemux, readablesize))

print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: amptk cluster -i %s -o out\n" % (catDemux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (catDemux)
