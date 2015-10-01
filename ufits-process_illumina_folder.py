#!/usr/bin/env python

#This script is a wrapper for -fastq_mergepairs from USEARCH8
import os, sys, argparse, shutil, subprocess, glob, math, logging
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

parser=argparse.ArgumentParser(prog='ufits-process_illumina_folder.py', usage="%(prog)s [options] -i folder",
    description='''Script that takes De-mulitplexed Illumina data from a folder and processes it for ufits (merge PE reads, strip primers, trim/pad to set length.''',
    epilog="""Written by Jon Palmer (2015) palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', dest='input', required=True, help='Folder of Illumina PE Data')
parser.add_argument('-o','--out', dest="out", default='ufits-results', help='Name for output folder')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='GTGARTCATCGAATCTTTG', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='TCCTCCGCTTATTGATATGC', help='Reverse Primer (ITS4)')
parser.add_argument('-n','--name_prefix', dest="prefix", default='R_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()      

FNULL = open(os.devnull, 'w')

usearch = args.usearch
try:
    subprocess.call([usearch, '--version'], stdout = FNULL, stderr = FNULL)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)

def convertSize(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Y', suffix) 

def myround(x, base=10):
    return int(base * round(float(x)/base))

def setupLogging(LOGNAME):
    global log
    stdoutformat = logging.Formatter(col.GRN+'%(asctime)s'+col.END+': %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
    log = logging.getLogger("myapp")
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)

#create directory and check for existing logfile
if not os.path.exists(args.out):
    os.makedirs(args.out)
    
log_name = os.path.join(args.out, 'ufits.log')
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print "-------------------------------------------------------"

#get filenames, store in list
#Illumina file names look like the following
#<sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz
filenames = []
for file in os.listdir(args.input):
    if file.endswith(".fastq"):
        filenames.append(file)

if len(filenames) % 2 != 0:
    print "Check your input files, they do not seem to be properly paired"
    os.exit(1)

uniq_names = []
fastq_for = []
fastq_rev = []
map = os.path.join(args.out, 'ufits-filenames.txt')
map_file = open(map, 'wb')
map_file.write("Name\t[i5]\t[i7]\tLane\tSet_num\n")
for item in sorted(filenames):
    if '_R1' in item:
        fastq_for.append(item)
    if '_R2' in item:
        fastq_rev.append(item)
    column = item.split("_")
    if column[0] not in uniq_names:
        uniq_names.append(column[0])
        if "-" in column[1]:
            barcode = column[1].split("-")
            try:
                map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], barcode[0], barcode[1], column[2], column[4].split(".",1)[0]))
            except IndexError:
                log.debug("Non-standard names detected, skipping mapping file")
                
        else:
            try:
                map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], column[1], "None", column[2], column[4].split(".",1)[0]))
            except IndexError:
                log.debug("Non-standard names detected, skipping mapping file")
map_file.close()

#loop through each set
for i in range(len(fastq_for)):
    name = fastq_for[i].split("_")[0]
    for_reads = os.path.join(args.input, fastq_for[i])
    rev_reads = os.path.join(args.input, fastq_rev[i])
    log.info("Working on reads from sample %s" % name)
    #get read length
    fp = open(for_reads)
    for i, line in enumerate(fp):
        if i == 1:
            read_length = len(line)
            read_length = myround(read_length)
        elif i > 2:
            break
    fp.close()

    #now trim the last bp off of the Illumina data (there for phasing, i.e. 250 bp reads are 251 bp)
    pretrim_R1 = os.path.join(args.out, 'pretrim_R1.fq')
    pretrim_R2 = os.path.join(args.out, 'pretrim_R2.fq')
    log.info("Merging Overlaping Pairs using USEARCH8")
    log.debug("%s -fastq_filter %s -fastq_trunclen %s -fastqout %s" % (usearch, for_reads, str(read_length), pretrim_R1))
    log.debug("%s -fastq_filter %s -fastq_trunclen %s -fastqout %s" % (usearch, rev_reads, str(read_length), pretrim_R2))
    subprocess.call([usearch, '-fastq_filter', for_reads, '-fastq_trunclen', str(read_length), '-fastqout', pretrim_R1], stdout = FNULL, stderr = FNULL)
    subprocess.call([usearch, '-fastq_filter', rev_reads, '-fastq_trunclen', str(read_length), '-fastqout', pretrim_R2], stdout = FNULL, stderr = FNULL)

    #next run USEARCH8 mergepe
    merge_out = os.path.join(args.out, 'merged.fq')
    skip_for = os.path.join(args.out, 'notmerged.R1.fq')
    log.debug("%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastqout_notmerged_fwd %s -fastq_truncqual 5 -fastq_allowmergestagger -minhsp 12" % (usearch, pretrim_R1, pretrim_R2, merge_out, skip_for))
    subprocess.call([usearch, '-fastq_mergepairs', for_reads, '-reverse', rev_reads, '-fastqout', merge_out, '-fastqout_notmerged_fwd', skip_for, '-fastq_truncqual', '5','-fastq_allowmergestagger','-minhsp', '12'], stdout = FNULL, stderr = FNULL)

    #now concatenate files for downstream pre-process_illumina.py script
    outname = name + '.fq'
    final_out = os.path.join(args.out, outname)
    out_file = open(final_out, 'wb')
    shutil.copyfileobj(open(merge_out,'rb'), out_file)
    shutil.copyfileobj(open(skip_for,'rb'), out_file)
    out_file.close()

    #clean and close up intermediate files
    os.remove(merge_out)
    os.remove(pretrim_R1)
    os.remove(pretrim_R2)
    os.remove(skip_for)

    log.info("Strip primers, trim/pad to %s bp\n" % args.trim_len)
    
    #now rest of script for demultiplexing here
    MAX_PRIMER_MISMATCHES = 2
    FileName = final_out
    FwdPrimer = args.F_primer
    RevPrimer = args.R_primer
    LabelPrefix = args.prefix
    SampleLabel = name
    MinLen = int(args.min_len)
    TrimLen = int(args.trim_len)
    demuxname = name + '.demux.fq'
    DemuxOut = os.path.join(args.out, demuxname)
    out_file = open(DemuxOut, 'w')
    
    SeqCount = 0
    OutCount = 0
    FwdPrimerMismatchCount = 0
    RevPrimerStrippedCount = 0
    TooShortCount = 0
    PadCount = 0

    PL = len(FwdPrimer)
    RevPrimer = revcomp_lib.RevComp(RevPrimer)
    
    print >> sys.stderr, "Foward primer: ", FwdPrimer
    print >> sys.stderr, "Rev comp'd rev primer: ", RevPrimer

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

        fastq.WriteRec(out_file, Label, Seq, Qual)

    fastq.ReadRecs(FileName, OnRec)
    progress.FileDone("%u reads, %u outupt, %u bad fwd primer, %u rev primer stripped, %u too short" % \
          (SeqCount, OutCount, FwdPrimerMismatchCount, RevPrimerStrippedCount, TooShortCount))

    print >> sys.stderr, "%10u seqs" % SeqCount
    print >> sys.stderr, "%10u fwd primer mismatches (%.1f%% discarded)" % (FwdPrimerMismatchCount, FwdPrimerMismatchCount*100.0/SeqCount)
    print >> sys.stderr, "%10u rev primer stripped (%.2f%% kept)" % (RevPrimerStrippedCount, RevPrimerStrippedCount*100.0/SeqCount)
    print >> sys.stderr, "%10u padded (%.2f%%)" % (PadCount, PadCount*100.0/SeqCount)
    print >> sys.stderr, "%10u too short (%.2f%%)" % (TooShortCount, TooShortCount*100.0/SeqCount)
    print >> sys.stderr, "%10u output (%.1f%%)\n" % (OutCount, OutCount*100.0/SeqCount)
    out_file.close()

print "-------------------------------------------------------"
#Now concatenate all of the demuxed files together
log.info("Concatenating Demuxed Files")

catDemux = os.path.join(args.out, 'ufits.demux.fq')
with open(catDemux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(args.out,'*.demux.fq')):
        if filename == catDemux:
            continue
        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)
            
log.info("Counting FASTQ Records")
num_lines = sum(1 for line in open(catDemux))
total = int(num_lines) / 4
log.info('{0:,}'.format(total) + ' reads processed')
log.info("Output file:  %s" % catDemux)

#get file size and issue warning if over 4.0 GB
filesize = os.path.getsize(catDemux)
readablesize = convertSize(filesize)
log.info("File size:  %s" % readablesize)
if filesize >= 4294967296:
    print col.WARN + "Warning, file is larger than 4 GB, you will need USEARCH 64 bit to cluster OTUs" + col.END
    
    