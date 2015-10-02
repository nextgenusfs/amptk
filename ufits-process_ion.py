#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, logging
from Bio import SeqIO
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

#get script path and barcode file name
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pgm_barcodes = os.path.join(script_path, 'lib', 'pgm_barcodes.fa')

parser=argparse.ArgumentParser(prog='ufits-process_ion.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Robert Edgar, modified by Jon Palmer (2015) palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', required=True, help='FASTQ file')
parser.add_argument('-o','--out', dest="out", default='ion', help='Base name for output')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='AGTGARTCATCGAATCTTTG', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='TCCTCCGCTTATTGATATGC', help='Reverse Primer (ITS4)')
parser.add_argument('-b','--list_barcodes', dest="barcodes", default='all', help='Enter Barcodes used separated by commas')
parser.add_argument('-n','--name_prefix', dest="prefix", default='R_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
parser.add_argument('--mult_samples', dest="multi", default='False', help='Combine multiple samples (i.e. FACE1)')
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

log_name = args.out + '.demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print "-------------------------------------------------------"

MAX_PRIMER_MISMATCHES = 2

FileName = args.fastq
FwdPrimer = args.F_primer
RevPrimer = args.R_primer
LabelPrefix = args.prefix
SampleLabel = args.multi
MinLen = int(args.min_len)
TrimLen = int(args.trim_len)

#get base name of input file
base = args.fastq.split(".")
base = base[0]

#get barcode list
barcode_file = base + ".barcodes_used.fa"
if os.path.isfile(barcode_file):
    os.remove(barcode_file)
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
BarcodeFileName = barcode_file
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

Barcodes = fasta.ReadSeqsDict(BarcodeFileName)

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

	OutCount += 1
	if SampleLabel == "False":
	    Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + BarcodeLabel + ";"
	else:
	    Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + SampleLabel + "_" + BarcodeLabel + ";"

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

#get file size and issue warning if over 4.0 GB
filesize = os.path.getsize(demuxname)
readablesize = convertSize(filesize)
log.info("File size:  %s" % readablesize)
print "-------------------------------------------------------"
if filesize >= 4294967296:
    if 'win32' in sys.platform:
        print "\nWarning, file is larger than 4 GB, you will need USEARCH 64 bit to cluster OTUs"
    else:
        print col.WARN + "\nWarning, file is larger than 4 GB, you will need USEARCH 64 bit to cluster OTUs" + col.END
else:
    if 'win32' in sys.platform:
        print "\nExample of next cmd: ufits cluster -i %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: BC_5)\n" % (demuxname)
    else:
        print col.WARN + "\nExample of next cmd: " + col.END + "ufits cluster -i %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: BC_5)\n" % (demuxname)

