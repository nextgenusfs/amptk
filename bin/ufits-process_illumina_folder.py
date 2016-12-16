#!/usr/bin/env python

#This script is a wrapper for -fastq_mergepairs from USEARCH8
import os, sys, argparse, shutil, subprocess, glob, math, logging, gzip, inspect, multiprocessing, itertools, re
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

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
parser.add_argument('-m','--mapping_file', help='Mapping file: QIIME format can have extra meta data columns')
parser.add_argument('--reads', dest="reads", default='paired', choices=['paired', 'forward'], help='PE or forward reads')
parser.add_argument('--read_length', default=300, type=int, help='Read length, i.e. 2 x 300 bp = 300')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer (ITS4)')
parser.add_argument('--require_primer', dest="primer", default='on', choices=['on', 'off'], help='Require Fwd primer to be present')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--rescue_forward', default='on', choices=['on', 'off'], help='Rescue Not-merged forward reads')
parser.add_argument('--min_len', default=50, type=int, help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default=250, type=int, help='Trim length for reads')
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH executable')
parser.add_argument('--cleanup', action='store_true', help='Delete all intermediate files')
args=parser.parse_args()      
    
def MergeReads(R1, R2, outname, read_length):
    usearch = args.usearch
    pretrim_R1 = outname + '.pretrim_R1.fq'
    pretrim_R2 = outname + '.pretrim_R2.fq'
    ufitslib.log.debug("Removing index 3prime bp 'A' from reads")    
    cmd = ['vsearch', '--fastq_filter', R1, '--fastq_trunclen', str(read_length), '--fastqout', pretrim_R1]
    ufitslib.runSubprocess(cmd, ufitslib.log)
    cmd = ['vsearch', '--fastq_filter', R2, '--fastq_trunclen', str(read_length), '--fastqout', pretrim_R2]
    ufitslib.runSubprocess(cmd, ufitslib.log)

    #next run USEARCH mergepe
    merge_out = outname + '.merged.fq'
    skip_for = outname + '.notmerged.R1.fq'
    ufitslib.log.debug("Now merging PE reads")
    cmd = [usearch, '-fastq_mergepairs', for_reads, '-reverse', rev_reads, '-fastqout', merge_out, '-fastqout_notmerged_fwd', skip_for,'-minhsp', '12','-fastq_maxdiffs', '8']
    ufitslib.runSubprocess(cmd, ufitslib.log)

    #now concatenate files for downstream pre-process_illumina.py script
    outname = outname + '.fq'
    final_out = os.path.join(args.out, outname)
    with open(final_out, 'w') as cat_file:
        shutil.copyfileobj(open(merge_out,'rU'), cat_file)
        if args.rescue_forward == 'on':
            shutil.copyfileobj(open(skip_for,'rU'), cat_file)
    
    #count output
    origcount = ufitslib.countfastq(R1)
    finalcount = ufitslib.countfastq(final_out)
    pct_out = finalcount / float(origcount)

    #clean and close up intermediate files
    os.remove(merge_out)
    os.remove(pretrim_R1)
    os.remove(pretrim_R2)
    os.remove(skip_for)
    return ufitslib.log.info('{0:,}'.format(finalcount) + ' reads passed ('+'{0:.1%}'.format(pct_out)+')')

def MatchesPrimer(Seq, Primer):
    return primer.MatchPrefix(Seq, Primer)

def processRead(input):
    #input is expected to be a FASTQ file
    #local variables that need to be previously declared: ForPrimer, RevPrimer
    Name = os.path.basename(input).split(".fq",-1)[0]
    DemuxOut = os.path.join(args.out, Name + '.demux.fq')
    counter = 1
    PL = len(FwdPrimer)
    with open(DemuxOut, 'w') as out:
        for title, seq, qual in FastqGeneralIterator(open(input)):
            #first thing is look for forward primer, if found trim it off
            Diffs = MatchesPrimer(seq, FwdPrimer)
            if args.primer == 'on':
                if Diffs > args.primer_mismatch:
                    continue
                else:
                    seq = seq[PL:]
                    qual = qual[PL:]
            elif args.primer == 'off':
                if Diffs <= args.primer_mismatch:
                    seq = seq[PL:]
                    qual = qual[PL:]
            #now look for reverse primer
            BestPosRev, BestDiffsRev = primer.BestMatch2(seq, RevPrimer, args.primer_mismatch)
            if BestPosRev > 0:  #reverse primer was found    
                #location to trim sequences, trim seqs
                Seq = seq[:BestPosRev]
                Qual = qual[:BestPosRev]
            #if full_length is passed, then only trim primers
            if not args.full_length:
                #got here if primers were found they were trimmed
                #now check seq length, pad if too short, trim if too long
                if len(Seq) < args.trim_len:
                    pad = args.trim_len - len(Seq)
                    Seq = Seq + pad*'N'
                    Qual = Qual +pad*'J'
                elif len(Seq) > args.trim_len:
                    Seq = Seq[:args.trim_len]
                    Qual = Qual[:args.trim_len]
            #got here, reads are primers trimmed and trim/padded, check length
            if len(Seq) < args.min_len:
                continue           
            #now fix header
            Title = 'R_'+str(counter)+';barcodelabel='+Name+';'
            #now write to file and bump counter
            counter += 1
            out.write("@%s\n%s\n+\n%s\n" % (Title, Seq, Qual))
                        
#sometimes people add slashes in the output directory, this could be bad, try to fix it
args.out = re.sub(r'\W+', '', args.out)
            
#create directory and check for existing logfile
if not os.path.exists(args.out):
    os.makedirs(args.out)
    
log_name = args.out+'.ufits-demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
ufitslib.SystemInfo()
#get version of ufits
usearch = args.usearch
ufitslib.versionDependencyChecks(usearch)

#check folder if files are gzipped, then gunzip them
#try to gunzip files
gzip_list = []
for file in os.listdir(args.input):
    if file.endswith(".fastq.gz"):
        gzip_list.append(file)
if gzip_list:
    ufitslib.log.info("Gzipped files detected, uncompressing")
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

#check for mapping file, if exists, then use names from first column only for filenames
if args.mapping_file:
    if not os.path.isfile(args.mapping_file):
        ufitslib.error("Mapping file is not valid: %s" % args.mapping_file)
        sys.exit(1)
    mapdata = ufitslib.parseMappingFileIllumina(args.mapping_file)
    #forward primer in first item in tuple, reverse in second
    sample_names = mapdata[0]
    FwdPrimer = mapdata[1]
    RevPrimer = mapdata[2]
    genericmapfile = args.mapping_file
    #loop through the files in the folder and get the ones in the sample_names lit
    filenames = []
    for file in os.listdir(args.input):
        if file.startswith(tuple(sample_names)):
            if file.endswith('.fastq'):
                filenames.append(file)
    
    if len(filenames) < 1:
        ufitslib.log.error("Found 0 valid files from mapping file. Mapping file SampleID must match start of filenames")
        sys.exit(1)

else: #if not then search through and find all the files you can in the folder
    '''get filenames, store in list, Illumina file names look like the following:
    <sample name>_<i5>-<i7>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz'''

    #now get the FASTQ files and proceed
    filenames = []
    for file in os.listdir(args.input):
        if file.endswith(".fastq"):
            filenames.append(file)
    #look up primer db otherwise default to entry
    if args.F_primer in ufitslib.primer_db:
        FwdPrimer = ufitslib.primer_db.get(args.F_primer)
    else:
        FwdPrimer = args.F_primer
    if args.R_primer in ufitslib.primer_db:
        RevPrimer = ufitslib.primer_db.get(args.R_primer)
    else:
        RevPrimer = args.R_primer

if len(filenames) % 2 != 0:
    print "Check your input files, they do not seem to be properly paired"
    sys.exit(1)

#check list for files, i.e. they need to have _R1 and _R2 in the filenames, otherwise throw exception
if not any('_R1' in x for x in filenames):
    ufitslib.log.error("Did not find valid FASTQ files.  Your files must have _R1 and _R2 in filename, rename your files and restart script.")
    sys.exit(1)

uniq_names = []
fastq_for = []
fastq_rev = []
sampleDict = {}
map = args.out + '.filenames.txt'
with open(map, 'w') as map_file:
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
                barcode = column[1].split("-")#looking here for the linker between i5 and i7 seqs
                i5 = barcode[0]
                i7 = barcode[1]
                try:
                    map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], i5, i7, column[2], column[4].split(".",1)[0]))
                except IndexError:
                    ufitslib.log.debug("Non-standard names detected, skipping mapping file")              
            else:
                i5 = column[1]
                i7 = "None"
                try:
                    map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], i5, i7, column[2], column[4].split(".",1)[0]))
                except IndexError:
                    ufitslib.log.debug("Non-standard names detected, skipping mapping file")
            if i7 != "None":
                sampleDict[column[0]] = i5+'-'+i7
            else:
                sampleDict[column[0]] = i5
#loop through each set and merge reads
if args.reads == 'paired':
    ufitslib.log.info("Merging Overlaping Pairs using USEARCH")

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
        if not file.endswith('.demux.fq'): #i don't want to demux the demuxed files.
            file = os.path.join(args.out, file)
            file_list.append(file)
if not args.full_length:
    ufitslib.log.info("Stripping primers and trim/pad to %s bp" % (args.trim_len))
else:
    ufitslib.log.info("Stripping primers and keeping only full length sequences")
ufitslib.log.info("splitting the job over %i cpus, but this may still take awhile" % (cpus))

#make sure primer is reverse complemented
RevPrimer = revcomp_lib.RevComp(RevPrimer)
#finally process reads over number of cpus
ufitslib.runMultiProgress(processRead, file_list, cpus)
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

if not args.mapping_file:
    #create a generic mappingfile for downstream processes
    genericmapfile = args.out + '.mapping_file.txt'
    ufitslib.CreateGenericMappingFileIllumina(sampleDict, FwdPrimer, revcomp_lib.RevComp(RevPrimer), genericmapfile)


#get file size
filesize = os.path.getsize(catDemux)
readablesize = ufitslib.convertSize(filesize)
ufitslib.log.info("Output file:  %s (%s)" % (catDemux, readablesize))
ufitslib.log.info("Mapping file: %s" % genericmapfile)
if args.cleanup:
    shutil.rmtree(args.out)
print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits cluster -i %s -o out\n" % (catDemux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "ufits cluster -i %s -o out\n" % (catDemux)
