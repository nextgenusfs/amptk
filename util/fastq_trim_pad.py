#!/usr/bin/env python

import sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.fasta as fasta
import lib.fastq as fastq
import lib.progress as progress
import lib.die as die

if len(sys.argv) < 2:
    print "Usage: fastq_trim_pad.py input.fastq trimlen"
    sys.exit()


FileName = sys.argv[1]
MinLen = 50
TrimLen = int(sys.argv[2])

SeqCount = 0
OutCount = 0
TooShortCount = 0
PadCount = 0


def OnRec(Label, Seq, Qual):
    global SeqCount, OutCount, TooShortCount, PadCount

    if SeqCount == 0:
        progress.InitFile(fastq.File)

    progress.File("%u reads, %u outupt, %u too short. %u padded" % \
      (SeqCount, OutCount, TooShortCount, PadCount))

    SeqCount += 1
    Seq = Seq
    Qual = Qual
    L = len(Seq)
    assert len(Qual) == L

    if L < MinLen:
        TooShortCount += 1
        return

    if L < TrimLen:
        PadCount += 1
        Seq = Seq + (TrimLen - L)*'N'
        Qual = Qual + (TrimLen - L)*'I'
        L = len(Seq)
        assert L == TrimLen
        assert len(Qual) == TrimLen

    if L > TrimLen:
        Seq = Seq[:TrimLen]
        Qual = Qual[:TrimLen]
        L = len(Seq)
    OutCount += 1
    
    assert L == TrimLen
    assert len(Qual) == TrimLen

    fastq.WriteRec(sys.stdout, Label, Seq, Qual)

fastq.ReadRecs(FileName, OnRec)
progress.FileDone("%u reads, %u outupt, %u too short" % \
      (SeqCount, OutCount, TooShortCount))

print >> sys.stderr, "%10u seqs" % SeqCount
print >> sys.stderr, "%10u padded (%.2f%%)" % (PadCount, PadCount*100.0/SeqCount)
print >> sys.stderr, "%10u too short (%.2f%%)" % (TooShortCount, TooShortCount*100.0/SeqCount)
print >> sys.stderr, "%10u output (%.1f%%)" % (OutCount, OutCount*100.0/SeqCount)
