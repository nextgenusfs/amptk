#!/usr/bin/env python

import sys, os, re, gzip, subprocess, argparse, inspect, logging, csv, shutil
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.fasta as fasta
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.amptklib as amptklib
from Bio.SeqIO.QualityIO import FastqGeneralIterator

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-fastq2sra.py', usage="%(prog)s [options] -i folder",
    description='''Script to split FASTQ file from Ion, 454, or Illumina by barcode sequence into separate files for submission to SRA.  This script can take the BioSample worksheet from NCBI and create an SRA metadata file for submission.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', dest='FASTQ', required=True, help='Input FASTQ file or folder')
parser.add_argument('-o','--out', dest='out', default="sra", help='Basename for output folder/files')
parser.add_argument('--min_len', default=50, type=int, help='Minimum length of read to keep')
parser.add_argument('-b','--barcode_fasta', dest='barcodes', help='Multi-fasta file containing barcodes used')
parser.add_argument('-s','--biosample', dest='biosample', help='BioSample file from NCBI')
parser.add_argument('-p','--platform', dest='platform', default='ion', choices=['ion', 'illumina', '454'], help='Sequencing platform')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer (ITS4)')
parser.add_argument('-n', '--names', help='CSV mapping file BC,NewName')
parser.add_argument('-d', '--description', help='Paragraph description for SRA metadata')
parser.add_argument('--force', action='store_true', help='Overwrite existing directory')
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


log_name = args.out + '.amptk-sra.log'
if os.path.isfile(log_name):
    os.remove(log_name)

amptklib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"

if args.platform != 'illumina' and not args.barcodes:
    amptklib.log.error("For ion, 454, or illumina2 datasets you must specificy a multi-fasta file containing barcodes with -b or --barcode_fasta")
    os._exit(1)

#initialize script, log system info
amptklib.log.info("Operating system: %s" % sys.platform)

#create output directory
if not os.path.exists(args.out):
    os.makedirs(args.out)
else:
    if not args.force:
        amptklib.log.error("Directory %s exists, add --force argument to overwrite" % args.out)
        os._exit(1)
    else:
        shutil.rmtree(args.out)
        os.makedirs(args.out)

if args.platform == 'illumina':
    #just need to get the correct .fastq.gz files into a folder by themselves
    #if illumina is selected, verify that args.fastq is a folder
    if not os.path.isdir(args.FASTQ):
        amptklib.log.error("%s is not a folder, for '--platform illumina', -i must be a folder containing raw reads" % (args.FASTQ))
        os._exit(1)
    rawlist = []
    filelist = []
    for file in os.listdir(args.FASTQ):
        if file.endswith(".fastq.gz"):
            rawlist.append(file)
    if len(rawlist) > 0:
        if not '_R2' in sorted(rawlist)[1]:
            amptklib.log.info("Found %i single files, copying to %s folder" % (len(rawlist), args.out))
            filelist = rawlist
            for file in rawlist:
                shutil.copyfile(os.path.join(args.FASTQ,file),(os.path.join(args.out,file)))
        else:
            amptklib.log.info("Found %i paired-end files, copying to %s folder" % (len(rawlist) / 2, args.out))
            for file in rawlist:
                shutil.copyfile(os.path.join(args.FASTQ,file),(os.path.join(args.out,file)))
                if '_R1' in file:
                    filelist.append(file)

else:
    #count FASTQ records in input
    amptklib.log.info("Loading FASTQ Records")
    total = amptklib.countfastq(args.FASTQ)
    size = amptklib.checkfastqsize(args.FASTQ)
    readablesize = amptklib.convertSize(size)
    amptklib.log.info('{0:,}'.format(total) + ' reads (' + readablesize + ')')

    #if --names given, load into dictonary
    if args.names:
        with open(args.names, 'rU') as input:
            reader = csv.reader(input)
            namesDict = {col[0]:col[1] for col in reader}
    else:
        amptklib.log.info("No names csv passed, using BC header names")

    #load barcode fasta file into dictonary
    Barcodes = {}
    files = []
    with open(args.barcodes, 'rU') as input:
        for line in input:
            if line.startswith('>'):
                if args.names:
                    name = namesDict.get(line[1:-1])
                    name = name + ".fastq"          
                else:
                    name = line[1:-1] + ".fastq"
                files.append(os.path.join(args.out,name))
                continue
            Barcodes[name]=line.strip()

    #ensure file names are unique        
    files = set(files)

    #this way will loop through the FASTQ file many times....not really what I want but it will work...
    runningTotal = 0
    with open(args.FASTQ, 'rU') as input:
        for title, seq, qual in FastqGeneralIterator(input):
            Barcode, BarcodeLabel = FindBarcode(seq)
            if Barcode == "": #if not found, move onto next record
                continue
            BarcodeLength = len(Barcode)
            seq = seq[BarcodeLength:]
            qual = qual[BarcodeLength:]
            #look for forward primer
            Diffs = MatchesPrimer(seq, FwdPrimer)
            if Diffs > 2:
                continue
            #if found, trim away primer
            seq = seq[len(FwdPrimer):]
            qual = qual[len(FwdPrimer):]
            #look for reverse primer, strip if found
            ReverseCompRev = revcomp_lib.RevComp(RevPrimer)
            BestPosRev, BestDiffsRev = primer.BestMatch2(seq, ReverseCompRev, 2)
            if BestPosRev > 0:
                seq = seq[:BestPosRev]
                qual = qual[:BestPosRev]
            #check size
            if len(seq) < args.min_len: #filter out sequences less than 50 bp.
                continue
            fileout = os.path.join(args.out, BarcodeLabel)
            with open(fileout, 'ab') as output:
                output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        Count = amptklib.countfastq(f)
        amptklib.log.info('{0:,}'.format(Count) + ' reads contained valid barcodes')
        runningTotal += Count

    amptklib.log.info('{0:,}'.format(runningTotal) + ' total reads for sra submission')



        
