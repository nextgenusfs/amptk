#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, logging, subprocess
from Bio import SeqIO
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.fasta as fasta
import lib.fastq as fastq
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.progress as progress
import lib.die as die

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-process_ion.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Robert Edgar, modified by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq','--sff', '--fasta', required=True, help='FASTQ/SFF/FASTA file')
parser.add_argument('-q','--qual', help='QUAL file (if -i is FASTA)')
parser.add_argument('-o','--out', dest="out", default='ion', help='Base name for output')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--barcode_fasta', default='pgm_barcodes.fa', help='FASTA file containing Barcodes (Names & Sequences)')
parser.add_argument('-b','--list_barcodes', dest="barcodes", default='all', help='Enter Barcodes used separated by commas')
parser.add_argument('-n','--name_prefix', dest="prefix", default='R_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
parser.add_argument('--mult_samples', dest="multi", default='False', help='Combine multiple samples (i.e. FACE1)')
parser.add_argument('--illumina', action='store_true', help='Input data is single file Illumina')
parser.add_argument('--ion', action='store_true', help='Input data is Ion Torrent')
parser.add_argument('--454', action='store_true', help='Input data is 454')
parser.add_argument('--reverse', help='Illumina reverse reads')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

def faqual2fastq(fasta, qual, fastq):
    global skipCount
    from Bio.SeqIO.QualityIO import PairedFastaQualIterator
    with open(fastq, 'w') as output:
        records = PairedFastaQualIterator(open(fasta), open(qual))
        for rec in records:
            try:
                SeqIO.write(rec, output, 'fastq')
            except ValueError:
                skipCount +1
    return skipCount
    
def convertSize(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Y', suffix) 

def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    else:
        stdoutformat = logging.Formatter(col.GRN+'%(asctime)s'+col.END+': %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)

log_name = args.out + '.demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print "-------------------------------------------------------"
#initialize script, log system info and usearch version
log.info("Operating system: %s" % sys.platform)

#if SFF file passed, convert to FASTQ with biopython
if args.fastq.endswith(".sff"):
    if args.barcode_fasta == 'pgm_barcodes.fa':
        log.error("You did not specify a --barcode_fasta, it is required for 454 data")
        os._exit(1)
    log.info("SFF input detected, converting to FASTQ")
    SeqIn = args.out + '.sff.extract.fastq'
    SeqIO.convert(args.fastq, "sff-trim", SeqIn, "fastq")
elif args.fastq.endswith(".fas") or args.fastq.endswith(".fasta") or args.fastq.endswith(".fa"):
    if not args.qual:
        log.error("FASTA input detected, however no QUAL file was given.  You must have FASTA + QUAL files")
        os._exit(1)
    else:
        if args.barcode_fasta == 'pgm_barcodes.fa':
            log.error("You did not specify a --barcode_fasta, it is required for 454 data")
            os._exit(1)
        SeqIn = args.out + '.fastq'
        log.info("FASTA + QUAL detected, converting to FASTQ")
        faqual2fastq(args.fastq, args.qual, SeqIn)
else:
    SeqIn = args.fastq

#check if illumina argument is passed, if so then run merge PE
if args.illumina:
    if args.barcode_fasta == 'pgm_barcodes.fa':
        log.error("You did not specify a --barcode_fasta, it is required this type of data")
        os._exit(1)
    if args.reverse:
        FNULL = open(os.devnull, 'w')
        #test for usearch
        usearch = args.usearch
        try:
            usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
        except OSError:
            log.warning("%s not found in your PATH, exiting." % usearch)
            os._exit(1)
        log.info("USEARCH version: %s" % usearch_test)
        
        #next run USEARCH8 mergepe
        SeqIn = args.out + '.merged.fq'
        log.info("Merging PE Illumina reads with USEARCH")
        log.debug("%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastq_truncqual 5 -fastq_maxdiffs 8 -minhsp 12" % (usearch, args.fastq, args.reverse, SeqIn))
        subprocess.call([usearch, '-fastq_mergepairs', args.fastq, '-reverse', args.reverse, '-fastqout', SeqIn, '-fastq_truncqual', '5','-minhsp', '12','-fastq_maxdiffs', '8'], stdout = FNULL, stderr = FNULL)
    else:
        log.info("Running UFITS on forward Illumina reads")
        SeqIn = args.fastq 

#look up primer db otherwise default to entry
if args.F_primer in ufitslib.primer_db:
    FwdPrimer = ufitslib.primer_db.get(args.F_primer)
else:
    FwdPrimer = args.F_primer
if args.R_primer in ufitslib.primer_db:
    RevPrimer = ufitslib.primer_db.get(args.R_primer)
else:
    RevPrimer = args.R_primer

#because we use an 'A' linker between barcode and primer sequence, add an A if ion is chemistry
if args.ion:
    FwdPrimer = 'A' + FwdPrimer

MAX_PRIMER_MISMATCHES = int(args.primer_mismatch)
FileName = SeqIn
LabelPrefix = args.prefix
SampleLabel = args.multi
MinLen = int(args.min_len)
TrimLen = int(args.trim_len)

#get base name of input file
base = args.fastq.split(".")
base = base[0]

#dealing with Barcodes
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

RevPrimer = revcomp_lib.RevComp(RevPrimer)
log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

SeqCount = 0
OutCount = 0
BarcodeMismatchCount = 0
FwdPrimerMismatchCount = 0
RevPrimerStrippedCount = 0
TooShortCount = 0
PadCount = 0
demuxname = args.out + '.demux.fq'
out_file = open(demuxname, 'w')

PL = len(FwdPrimer)

Barcodes = fasta.ReadSeqsDict(barcode_file)
BarcodeCount = {}
def MatchesPrimer(Seq, Primer):
    return primer.MatchPrefix(Seq, Primer)

def FindBarcode(Seq):
    global Barcodes
    for BarcodeLabel in Barcodes.keys():
        Barcode = Barcodes[BarcodeLabel]
        if Seq.startswith(Barcode):
            return Barcode, BarcodeLabel
    return "", ""

def OnRec(Label, Seq, Qual):
    global PL, LabelPrefix, Barcode, SeqCount, OutCount, TooShortCount, PadCount
    global BarcodeMismatchCount, FwdPrimerMismatchCount, RevPrimerStrippedCount
    global FwdPrimer, RevPrimer
    global BarcodeCount

    if SeqCount == 0:
        progress.InitFile(fastq.File)

    progress.File("%u reads, %u outupt, %u bad barcode, %u bad fwd primer, %u rev primer stripped, %u too short. %u padded" % \
      (SeqCount, OutCount, BarcodeMismatchCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount, PadCount))

    SeqCount += 1
    Barcode, BarcodeLabel = FindBarcode(Seq)
    if Barcode == "":
        BarcodeMismatchCount += 1
        return
          
    BarcodeLength = len(Barcode)
    Seq = Seq[BarcodeLength:]
    Qual = Qual[BarcodeLength:]

    Diffs = MatchesPrimer(Seq, FwdPrimer)
    if Diffs > MAX_PRIMER_MISMATCHES:
        FwdPrimerMismatchCount += 1
        return

    # Strip fwd primer
    Seq = Seq[PL:]
    Qual = Qual[PL:]

    BestPosRev, BestDiffsRev = primer.BestMatch2(Seq, RevPrimer, MAX_PRIMER_MISMATCHES)
    if BestPosRev > 0:
        # Strip rev primer
        RevPrimerStrippedCount += 1
        StrippedSeq = Seq[:BestPosRev]
        StrippedQual = Qual[:BestPosRev]

        # correctness checks
        if 1:
            Tail = Seq[BestPosRev:]
            Diffs2 = primer.MatchPrefix(Tail, RevPrimer)
            if Diffs2 != BestDiffsRev:
                print >> sys.stderr
                print >> sys.stderr, " Seq=" + Seq
                print >> sys.stderr, "Tail=" + Tail
                print >> sys.stderr, "RevP=" + RevPrimer
                die.Die("BestPosRev %u Diffs2 %u BestDiffsRev %u" % (BestPosRev, Diffs2, BestDiffsRev))
            assert StrippedSeq + Tail == Seq

        Seq = StrippedSeq
        Qual = StrippedQual

        L = len(Seq)
        assert len(Qual) == L

        if L < MinLen:
            return

        if L < TrimLen:
            PadCount += 1
            Seq = Seq + (TrimLen - L)*'N'
            Qual = Qual + (TrimLen - L)*'I'
            L = len(Seq)
            assert L == TrimLen
            assert len(Qual) == TrimLen

    L = len(Seq)
    if L < TrimLen:
        TooShortCount += 1
        return

    if L > TrimLen:
        Seq = Seq[:TrimLen]
        Qual = Qual[:TrimLen]
        L = len(Seq)
        
    OutCount += 1
    
    if BarcodeLabel not in BarcodeCount:
        BarcodeCount[BarcodeLabel] = 1
    else:
        BarcodeCount[BarcodeLabel] += 1
    
    if SampleLabel == "False":
        Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + BarcodeLabel + ";"
    else:
        Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + SampleLabel + "_" + BarcodeLabel + ";"
    
    assert L == TrimLen
    assert len(Qual) == TrimLen

    fastq.WriteRec(out_file, Label, Seq, Qual)

fastq.ReadRecs(FileName, OnRec)
progress.FileDone("%u reads, %u outupt, %u bad barcode, %u bad fwd primer, %u rev primer stripped, %u too short" % \
      (SeqCount, OutCount, BarcodeMismatchCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount))

log.info("Stats for demuxing: \
\n%10u seqs \
\n%10u barcode mismatches \
\n%10u fwd primer mismatches (%.1f%% discarded) \
\n%10u rev primer stripped (%.2f%% kept) \
\n%10u padded (%.2f%%) \
\n%10u too short (%.2f%%) \
\n%10u output (%.1f%%)" % \
(SeqCount, BarcodeMismatchCount, FwdPrimerMismatchCount, FwdPrimerMismatchCount*100.0/SeqCount, RevPrimerStrippedCount, RevPrimerStrippedCount*100.0/SeqCount, PadCount, PadCount*100.0/SeqCount, TooShortCount, TooShortCount*100.0/SeqCount, OutCount, OutCount*100.0/SeqCount))
out_file.close()

#now let's count the barcodes found and count the number of times they are found.
barcode_counts = "%10s:  %s" % ('Sample', 'Count')
for key,value in natsorted(BarcodeCount.iteritems()):
    barcode_counts += "\n%10s:  %s" % (key, str(value))
log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

#get file size
filesize = os.path.getsize(demuxname)
readablesize = convertSize(filesize)
log.info("Output file: %s (%s)" % (demuxname, readablesize))

print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits cluster -i %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: BC_5)\n" % (demuxname)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "ufits cluster -i %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: BC_5)\n" % (demuxname)

