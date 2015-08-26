#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import lib.fasta as fasta
import lib.fastq as fastq
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.progress as progress
import lib.die as die

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

parser=argparse.ArgumentParser(prog='ufits-process_illumina.py', usage="%(prog)s [options] file.fastq > out.fastq\n%(prog)s -h for help menu",
    description='''Script strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Robert Edgar, modified slightly by Jon Palmer (2015) palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('fastq', help='FASTQ file')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='GTGARTCATCGAATCTTTG', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='TCCTCCGCTTATTGATATGC', help='Reverse Primer (ITS4)')
parser.add_argument('-n','--name_prefix', dest="prefix", default='Reads_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
parser.add_argument('-b','--barcode', dest="barcode", default="filename", help='Name of Sample for downstream processing')
args=parser.parse_args()

MAX_PRIMER_MISMATCHES = 2

FileName = args.fastq
FwdPrimer = args.F_primer
RevPrimer = args.R_primer
LabelPrefix = args.prefix
if args.barcode == "filename":
    SampleLabel = FileName.split(".")[0]
else:
    SampleLabel = args.barcode
MinLen = int(args.min_len)
TrimLen = int(args.trim_len)

print >> sys.stderr, "Foward primer: ", FwdPrimer
print >> sys.stderr, "Rev comp'd rev primer: ", RevPrimer

SeqCount = 0
OutCount = 0
FwdPrimerMismatchCount = 0
RevPrimerStrippedCount = 0
TooShortCount = 0
PadCount = 0

PL = len(FwdPrimer)
RevPrimer = revcomp_lib.RevComp(RevPrimer)

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

    assert L == TrimLen
    assert len(Qual) == TrimLen

    fastq.WriteRec(sys.stdout, Label, Seq, Qual)

fastq.ReadRecs(FileName, OnRec)
progress.FileDone("%u reads, %u outupt, %u bad fwd primer, %u rev primer stripped, %u too short" % \
      (SeqCount, OutCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount))

print >> sys.stderr, "%10u seqs" % SeqCount
print >> sys.stderr, "%10u fwd primer mismatches (%.1f%% discarded)" % (FwdPrimerMismatchCount, FwdPrimerMismatchCount*100.0/SeqCount)
print >> sys.stderr, "%10u rev primer stripped (%.2f%% kept)" % (RevPrimerStrippedCount, RevPrimerStrippedCount*100.0/SeqCount)
print >> sys.stderr, "%10u padded (%.2f%%)" % (PadCount, PadCount*100.0/SeqCount)
print >> sys.stderr, "%10u too short (%.2f%%)" % (TooShortCount, TooShortCount*100.0/SeqCount)
print >> sys.stderr, "%10u output (%.1f%%)" % (OutCount, OutCount*100.0/SeqCount)
