#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import object
import sys, os, argparse, shutil, logging, inspect
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
class col(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-process_illumina.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Robert Edgar, modified by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', required=True, help='FASTQ file')
parser.add_argument('-o','--out', dest="out", default='filename', help='Base name for output')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='GTGARTCATCGAATCTTTG', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='TCCTCCGCTTATTGATATGC', help='Reverse Primer (ITS4)')
parser.add_argument('-n','--name_prefix', dest="prefix", default='Reads_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
args=parser.parse_args()

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

FileName = args.fastq
FwdPrimer = args.F_primer
RevPrimer = args.R_primer
LabelPrefix = args.prefix
if args.out == "filename":
    SampleLabel = FileName.split(".")[0]
else:
    SampleLabel = args.out
MinLen = int(args.min_len)
TrimLen = int(args.trim_len)

log_name = SampleLabel + '.demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print("-------------------------------------------------------")


MAX_PRIMER_MISMATCHES = 2

demuxname = SampleLabel + '.demux.fq'
out_file = open(demuxname, 'w')

SeqCount = 0
OutCount = 0
FwdPrimerMismatchCount = 0
RevPrimerStrippedCount = 0
TooShortCount = 0
PadCount = 0

PL = len(FwdPrimer)
RevPrimer = revcomp_lib.RevComp(RevPrimer)

log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

def MatchesPrimer(Seq, Primer):
    return primer.MatchPrefix(Seq, Primer)

def OnRec(Label, Seq, Qual):
    global PL, LabelPrefix, SeqCount, OutCount, TooShortCount, PadCount
    global FwdPrimerMismatchCount, RevPrimerStrippedCount
    global FwdPrimer, RevPrimer

    if SeqCount == 0:
        progress.InitFile(fastq.File)

    progress.File("%u reads, %u outupt, %u bad fwd primer, %u rev primer stripped, %u too short. %u padded" % \
      (SeqCount, OutCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount, PadCount))

    SeqCount += 1
    Seq = Seq
    Qual = Qual
    Diffs = MatchesPrimer(Seq, FwdPrimer)
    if Diffs > MAX_PRIMER_MISMATCHES:
        FwdPrimerMismatchCount += 1
        return

    OutCount += 1
    Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + SampleLabel + ";"

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
                print(file=sys.stderr)
                print(" Seq=" + Seq, file=sys.stderr)
                print("Tail=" + Tail, file=sys.stderr)
                print("RevP=" + RevPrimer, file=sys.stderr)
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

    assert L == TrimLen
    assert len(Qual) == TrimLen

    fastq.WriteRec(out_file, Label, Seq, Qual)

fastq.ReadRecs(FileName, OnRec)
progress.FileDone("%u reads, %u outupt, %u bad fwd primer, %u rev primer stripped, %u too short" % \
      (SeqCount, OutCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount))

log.info("Stats for demuxing: \
\n%10u seqs \
\n%10u fwd primer mismatches (%.1f%% discarded) \
\n%10u rev primer stripped (%.2f%% kept) \
\n%10u padded (%.2f%%) \
\n%10u too short (%.2f%%) \
\n%10u output (%.1f%%)" % \
(SeqCount, FwdPrimerMismatchCount, FwdPrimerMismatchCount*100.0/SeqCount, RevPrimerStrippedCount, RevPrimerStrippedCount*100.0/SeqCount, PadCount, PadCount*100.0/SeqCount, TooShortCount, TooShortCount*100.0/SeqCount, OutCount, OutCount*100.0/SeqCount))
out_file.close()

#get file size and issue warning if over 4.0 GB
filesize = os.path.getsize(demuxname)
readablesize = convertSize(filesize)
log.info("File size:  %s" % readablesize)
print("-------------------------------------------------------")
if filesize >= 4294967296:
    print("\nWarning, file is larger than 4 GB, you will need USEARCH 64 bit to cluster OTUs")
else:
    print("\nExample of next cmd: ufits-OTU_cluster.py -f %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: BC_5)\n" % (demuxname))
