#!/usr/bin/env python

import sys, os, argparse, inspect, logging, shutil
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.fasta as fasta
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.amptklib as amptklib
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import edlib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-split_samples_fastq.py', usage="%(prog)s [options] -i file.fastq",
    description='''Script to split FASTQ into indifidual files.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', dest='FASTQ', required=True, help='Input FASTQ file or folder')
parser.add_argument('-o','--out', dest='out', default="amptk-split", help='Basename for output folder/files')
parser.add_argument('--min_len', default=50, type=int, help='Minimum length of read to keep')
parser.add_argument('-b','--barcode_fasta', dest='barcodes', required=True, help='Multi-fasta file containing barcodes used')
parser.add_argument('-p','--platform', dest='platform', default='ion', choices=['ion', 'illumina', '454'], help='Sequencing platform')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer (ITS4)')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--force', action='store_true', help='overwrite output folder')
args=parser.parse_args()


#check if name is in primer_db, else use input value
if args.F_primer in amptklib.primer_db:
    FwdPrimer = amptklib.primer_db.get(args.F_primer)
else:
    FwdPrimer = args.F_primer

#check if name is in primer_db, else use input value
if args.R_primer in amptklib.primer_db:
    RevPrimer = amptklib.primer_db.get(args.R_primer)
else:
    RevPrimer = args.R_primer
    
ReverseCompRev = revcomp_lib.RevComp(RevPrimer)
#add the linker for ion primers
if args.platform == 'ion':
    FwdPrimer = 'A' + FwdPrimer

def FindBarcode(Seq):
    global Barcodes
    for BarcodeLabel in Barcodes.keys():
        Barcode = Barcodes[BarcodeLabel]
        if Seq.startswith(Barcode):
            return Barcode, BarcodeLabel
    return "", ""

def MatchesPrimer(Seq, Primer):
    return primer.MatchPrefix(Seq, Primer)

log_name = args.out + '.amptk-split.log'
if os.path.isfile(log_name):
    os.remove(log_name)

amptklib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info
amptklib.log.info("Operating system: %s, %s" % (sys.platform, amptklib.get_version()))

#create output directory
if not os.path.exists(args.out):
    os.makedirs(args.out)
else:
    if not args.force:
        amptklib.log.error("Directory %s exists, add --force argument to overwrite" % args.out)
    else:
        shutil.rmtree(args.out)
        os.makedirs(args.out)

#count FASTQ records in input
amptklib.log.info("Loading FASTQ Records")
total = amptklib.countfastq(args.FASTQ)
size = amptklib.checkfastqsize(args.FASTQ)
readablesize = amptklib.convertSize(size)
amptklib.log.info('{0:,}'.format(total) + ' reads (' + readablesize + ')')

#load barcode fasta file into dictonary
Barcodes = {}
with open(args.barcodes, 'rU') as input:
    for line in input:
        if line.startswith('>'):
            name = line[1:-1] + ".fastq"
            continue
        Barcodes[name]=line.strip()

amptklib.log.info("Looking for %i barcodes and trimming primers\nFwdPrimer: %s\nRevPrimer: %s" % (len(Barcodes), FwdPrimer, RevPrimer))

#this will loop through FASTQ file once, splitting those where barcodes are found, and primers trimmed
runningTotal = 0
badPrimer = 0
noBC = 0
BC = 0
trim = len(FwdPrimer)
with open(args.FASTQ, 'rU') as input:
    for title, seq, qual in FastqGeneralIterator(input):
        Barcode, BarcodeLabel = FindBarcode(seq)
        if Barcode == "": #if not found, move onto next record
            noBC += 1
            continue
        BC += 1
        BarcodeLength = len(Barcode)
        seq = seq[BarcodeLength:]
        qual = qual[BarcodeLength:]
        #look for forward primer
        Diffs = MatchesPrimer(seq, FwdPrimer)
        if Diffs > args.primer_mismatch:
            badPrimer += 1
            continue
        #if found, trim away primer
        seq = seq[trim:]
        qual = qual[trim:]
        #look for reverse primer, strip if found
        BestPosRev, BestDiffsRev = primer.BestMatch2(seq, ReverseCompRev, args.primer_mismatch)
        if BestPosRev > 0:
            seq = seq[:BestPosRev]
            qual = qual[:BestPosRev]
        #check size
        if len(seq) < args.min_len: #filter out sequences less than 50 bp.
            continue
        runningTotal += 1
        fileout = os.path.join(args.out, BarcodeLabel)
        with open(fileout, 'ab') as output:
            output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))    
amptklib.log.info('{0:,}'.format(runningTotal) + ' total reads with valid barcode')
amptklib.log.info('{0:,}'.format(BC) + ' barcode found')
amptklib.log.info('{0:,}'.format(noBC) + ' no barcode found')
amptklib.log.info('{0:,}'.format(badPrimer) + ' no fwd primer found')