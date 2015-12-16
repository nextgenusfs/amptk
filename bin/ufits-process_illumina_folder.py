#!/usr/bin/env python

#This script is a wrapper for -fastq_mergepairs from USEARCH8
import os, sys, argparse, shutil, subprocess, glob, math, logging, gzip, inspect, multiprocessing, itertools
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
from Bio import SeqIO

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-process_illumina_folder.py', usage="%(prog)s [options] -i folder",
    description='''Script that takes De-mulitplexed Illumina data from a folder and processes it for ufits (merge PE reads, strip primers, trim/pad to set length.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', dest='input', required=True, help='Folder of Illumina Data')
parser.add_argument('-o','--out', dest="out", default='ufits-data', help='Name for output folder')
parser.add_argument('--reads', dest="reads", default='paired', choices=['paired', 'forward'], help='PE or forward reads')
parser.add_argument('--read_length', default=300, type=int, help='Read length, i.e. 2 x 300 bp = 300')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer (ITS4)')
parser.add_argument('--require_primer', dest="primer", default='on', choices=['on', 'off'], help='Require Fwd primer to be present')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--rescue_forward', action="store_true", help='Rescue Not-merged forward reads')
parser.add_argument('-n','--name_prefix', dest="prefix", default='R_', help='Prefix for renaming reads')
parser.add_argument('-m','--min_len', default='50', help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default='250', help='Trim length for reads')
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()      

#look up primer db otherwise default to entry
if args.F_primer in ufitslib.primer_db:
    FwdPrimer = ufitslib.primer_db.get(args.F_primer)
else:
    FwdPrimer = args.F_primer
if args.R_primer in ufitslib.primer_db:
    RevPrimer = ufitslib.primer_db.get(args.R_primer)
else:
    RevPrimer = args.R_primer
    
def MergeReads(R1, R2, outname, read_length):
    usearch = args.usearch
    pretrim_R1 = outname + '.pretrim_R1.fq'
    pretrim_R2 = outname + '.pretrim_R2.fq'
    ufitslib.log.debug("%s -fastq_filter %s -fastq_trunclen %s -fastqout %s" % (usearch, R1, str(read_length), pretrim_R1))
    ufitslib.log.debug("%s -fastq_filter %s -fastq_trunclen %s -fastqout %s" % (usearch, R2, str(read_length), pretrim_R2))
    subprocess.call([usearch, '-fastq_filter', R1, '-fastq_trunclen', str(read_length), '-fastqout', pretrim_R1], stdout = FNULL, stderr = FNULL)
    subprocess.call([usearch, '-fastq_filter', R2, '-fastq_trunclen', str(read_length), '-fastqout', pretrim_R2], stdout = FNULL, stderr = FNULL)

    #next run USEARCH8 mergepe
    merge_out = outname + '.merged.fq'
    skip_for = outname + '.notmerged.R1.fq'
    ufitslib.log.debug("%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastqout_notmerged_fwd %s -fastq_truncqual 5 -fastq_maxdiffs 8 -minhsp 12" % (usearch, pretrim_R1, pretrim_R2, merge_out, skip_for))
    subprocess.call([usearch, '-fastq_mergepairs', for_reads, '-reverse', rev_reads, '-fastqout', merge_out, '-fastqout_notmerged_fwd', skip_for, '-fastq_truncqual', '5','-minhsp', '12','-fastq_maxdiffs', '8'], stdout = FNULL, stderr = FNULL)

    #now concatenate files for downstream pre-process_illumina.py script
    outname = outname + '.fq'
    final_out = os.path.join(args.out, outname)
    cat_file = open(final_out, 'w')
    shutil.copyfileobj(open(merge_out,'rU'), cat_file)
    if args.rescue_forward:
        shutil.copyfileobj(open(skip_for,'rU'), cat_file)
    cat_file.close()

    #clean and close up intermediate files
    os.remove(merge_out)
    os.remove(pretrim_R1)
    os.remove(pretrim_R2)
    os.remove(skip_for)

def MatchesPrimer(Seq, Primer):
    return primer.MatchPrefix(Seq, Primer)
    
def ProcessReads(records):
    OutCount = 0
    MAX_PRIMER_MISMATCHES = int(args.primer_mismatch)
    LabelPrefix = args.prefix
    MinLen = int(args.min_len)
    TrimLen = int(args.trim_len)
    PL = len(FwdPrimer)
    revPrimer = revcomp_lib.RevComp(RevPrimer)
    for rec in records:
        OutCount += 1
        rec.id = LabelPrefix + str(OutCount) + ";barcodelabel=" + name + ";"
        rec.name = ""
        rec.description = ""
        #turn sequence into string for matching
        Seq = str(rec.seq)
        Diffs = MatchesPrimer(Seq, FwdPrimer)
        if args.primer == "on":
            if Diffs > MAX_PRIMER_MISMATCHES:
                continue
            # Strip fwd primer from rec
            rec = rec[PL:]
        elif args.primer == "off":
            if Diffs < MAX_PRIMER_MISMATCHES:
                # Strip fwd primer from rec
                rec = rec[PL:]
                
        #turn seq into str again
        Seq = str(rec.seq)
        #look for reverse primer
        BestPosRev, BestDiffsRev = primer.BestMatch2(Seq, revPrimer, MAX_PRIMER_MISMATCHES)
        if BestPosRev > 0:
            # Strip rev primer from rec.seq
            rec = rec[:BestPosRev]
            #check length       
            L = len(rec.seq)
            if L < MinLen:
                continue
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
            #check length
            L = len(rec.seq)
            #truncate down to trim length
            if L >= TrimLen:
                rec = rec[:TrimLen]        
                yield rec

def worker(file):
    global name
    name = file.split(".",-1)[0]
    name = name.split("/",1)[1]
    demuxname = name + '.demux.fq'
    DemuxOut = os.path.join(args.out, demuxname)
    with open(DemuxOut, 'w') as output:
        with open(file, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'fastq')
            SeqIO.write(ProcessReads(SeqRecords), output, 'fastq')
            
#create directory and check for existing logfile
if not os.path.exists(args.out):
    os.makedirs(args.out)
    
log_name = os.path.join(args.out, 'ufits.log')
if os.path.isfile(log_name):
    os.remove(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
ufitslib.log.info("Operating system: %s" % sys.platform)
usearch = args.usearch
try:
    usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
except OSError:
    ufitslib.log.warning("%s not found in your PATH, exiting." % usearch)
    os._exit(1)
ufitslib.log.info("USEARCH version: %s" % usearch_test)

'''get filenames, store in list, Illumina file names look like the following:
<sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz'''
#try to gunzip files
gzip_list = []
for file in os.listdir(args.input):
    if file.endswith(".fastq.gz"):
        gzip_list.append(file)
if gzip_list:
    ufitslib.log.info("Gzipped files detected, uncompressing")

#check list for valid filenames they need to have _R1 and _R2, otherwise through exception
if gzip_list and '_R1' not in gzip_list[0]:
    ufitslib.log.error("Did not find valid FASTQ files.  Your files must have _R1 and _R2 in filename, rename your files and restart script.")
    os._exit(1)

if gzip_list:
    for file in gzip_list:
        ufitslib.log.debug("Uncompressing %s" % file)
        OutName = os.path.join(args.input, os.path.splitext(file)[0])
        InFile = gzip.open(os.path.join(args.input, file), 'rU')
        ReadFile = InFile.read()
        OutFile = open(OutName, 'w')
        OutFile.write(ReadFile)
        OutFile.close()
        InFile.close()
        os.remove(os.path.join(args.input, file)) #remove .gz file    

#now get the FASTQ files and proceed
filenames = []
for file in os.listdir(args.input):
    if file.endswith(".fastq"):
        filenames.append(file)

if len(filenames) % 2 != 0:
    print "Check your input files, they do not seem to be properly paired"
    os._exit(1)

#check list for files, i.e. they need to have _R1 and _R2 in the filenames, otherwise through exception
if '_R1' not in filenames[0]:
    ufitslib.log.error("Did not find valid FASTQ files.  Your files must have _R1 and _R2 in filename, rename your files and restart script.")
    os._exit(1)

uniq_names = []
fastq_for = []
fastq_rev = []
map = args.out + '-filenames.txt'
map_file = open(map, 'w')
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
                ufitslib.log.debug("Non-standard names detected, skipping mapping file")
                
        else:
            try:
                map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], column[1], "None", column[2], column[4].split(".",1)[0]))
            except IndexError:
                ufitslib.log.debug("Non-standard names detected, skipping mapping file")
map_file.close()

#loop through each set and merge reads
if args.reads == 'paired':
    ufitslib.log.info("Merging Overlaping Pairs using USEARCH8")
for i in range(len(fastq_for)):
    name = fastq_for[i].split("_")[0]
    outname = name + '.fq'
    if os.path.isfile(os.path.join(args.out, outname)):
        ufitslib.log.info("Output for %s detected, skipping files" % outname)
        continue
    for_reads = os.path.join(args.input, fastq_for[i])
    rev_reads = os.path.join(args.input, fastq_rev[i])
    ufitslib.log.info("working on sample %s" % name)
    if args.reads == 'paired':
        #get read length
        fp = open(for_reads)
        for i, line in enumerate(fp):
            if i == 1:
                read_length = len(line)
                read_length = ufitslib.myround(read_length)
            elif i > 2:
                break
        fp.close()

        if args.read_length != read_length:
            ufitslib.log.warning("Measured read length (%i bp) does not equal %i bp, proceeding with larger value" % (read_length, args.read_length))
            if args.read_length > read_length:
                read_length = args.read_length
            else:
                read_length = read_length
                
        MergeReads(for_reads, rev_reads, name, read_length)
    else:
        shutil.copy(for_reads, os.path.join(args.out, outname))
    
#Now all the data is in folder args.out that needs to be de-multiplexed
if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

#get list of files to demux
file_list = []
for file in os.listdir(args.out):
    if file.endswith(".fq"):
        file = os.path.join(args.out, file)
        file_list.append(file)
ufitslib.log.info("Stripping primers and trim/pad to %s bp" % (args.trim_len))
ufitslib.log.info("splitting the job over %i cpus, but this may still take awhile" % (cpus))

#parallize over number of cpus
p = multiprocessing.Pool(cpus)
for f in file_list:
    p.apply_async(worker, [f])
p.close()
p.join()

print "-------------------------------------------------------"
#Now concatenate all of the demuxed files together
ufitslib.log.info("Concatenating Demuxed Files")

catDemux = args.out + '.demux.fq'
with open(catDemux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(args.out,'*.demux.fq')):
        if filename == catDemux:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)
            
ufitslib.log.info("Counting FASTQ Records")
total = ufitslib.countfastq(catDemux)
ufitslib.log.info('{0:,}'.format(total) + ' reads processed')

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
ufitslib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

#get file size
filesize = os.path.getsize(catDemux)
readablesize = ufitslib.convertSize(filesize)
ufitslib.log.info("Output file:  %s (%s)" % (catDemux, readablesize))
print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits cluster -i %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: spike)\n" % (catDemux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "ufits cluster -i %s -o out --uchime_ref ITS2 --mock <mock BC name> (test data: spike)\n" % (catDemux)
