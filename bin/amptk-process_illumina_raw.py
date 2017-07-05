#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, logging, subprocess, multiprocessing, glob, itertools, re, edlib
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.fasta as fasta
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.amptklib as amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-process_illumina_raw.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--forward', dest='fastq', required=True, help='Illumina FASTQ R1 reads')
parser.add_argument('-r', '--reverse', required=True, help='Illumina FASTQ R2 reads')
parser.add_argument('-i', '--index', required=True, help='Illumina FASTQ index reads')
parser.add_argument('-m', '--mapping_file', required=True, help='QIIME-like mapping tool')
parser.add_argument('--read_length', type=int, help='Read length, i.e. 2 x 300 bp = 300')
parser.add_argument('-o','--out', dest="out", default='illumina_out', help='Base name for output')
parser.add_argument('--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer')
parser.add_argument('--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--barcode_mismatch', default=2, type=int, help='Number of mis-matches in barcode')
parser.add_argument('--barcode_fasta', default='pgm_barcodes.fa', help='FASTA file containing Barcodes (Names & Sequences)')
parser.add_argument('--require_primer', dest="primer", default='off', choices=['on', 'off'], help='Require Fwd primer to be present')
parser.add_argument('--rescue_forward', default='on', choices=['on', 'off'], help='Rescue Not-merged forward reads')
parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default=300, type=int, help='Trim length for reads')
parser.add_argument('-p','--pad', default='off', choices=['on', 'off'], help='Pad with Ns to a set length')
parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
args=parser.parse_args()

def processRead(input):
    #input is expected to be a FASTQ file
    #local variables that need to be previously declared: ForPrimer, RevPrimer, indexReads, discard
    Name = os.path.basename(input).split(".fq",-1)[0]
    DemuxOut = os.path.join(tmpdir, Name + '.demux.fq')
    Sample = Name.split('_')[0]
    StatsOut = os.path.join(tmpdir, Name+'.stats')
    Total = 0
    NoBC = 0
    NoPrimer = 0
    TooShort = 0
    RevPrimerFound = 0
    ValidSeqs = 0
    with open(StatsOut, 'w') as counts:
        with open(DemuxOut, 'w') as out:
            for title, seq, qual in FastqGeneralIterator(open(input)):
                Total += 1
                #check if in discard
                readID = title.split(' ')[0]
                if readID in discard:
                    NoBC += 1
                    continue
                #first thing is look for forward primer, if found trim it off
                foralign = edlib.align(FwdPrimer, seq, mode="HW", k=args.primer_mismatch)
                #if require primer is on make finding primer in amplicon required if amplicon is larger than read length
                #if less than read length, can't enforce primer because could have been trimmed via staggered trim in fastq_mergepairs
                if args.primer == 'on' and len(seq) > ReadLen:
                    if foralign["editDistance"] < 0:
                        NoPrimer += 1
                        continue
                    ForCutPos = foralign["locations"][0][1]+1
                    Seq = seq[ForCutPos:]
                    Qual = qual[ForCutPos:]
                else:
                    if foralign["editDistance"] >= 0:
                        ForCutPos = foralign["locations"][0][1]+1
                        Seq = seq[ForCutPos:]
                        Qual = qual[ForCutPos:]
                    else:
                        NoPrimer += 1
                        Seq = seq
                        Qual = qual
                #now look for reverse primer
                revalign = edlib.align(RevPrimer, Seq, mode="HW", task="locations", k=args.primer_mismatch)
                if revalign["editDistance"] >= 0:
                    RevPrimerFound += 1
                    RevCutPos = revalign["locations"][0][0]
                    #location to trim sequences, trim seqs
                    Seq = Seq[:RevCutPos]
                    Qual = Qual[:RevCutPos]
                else:
                    if args.full_length and len(Seq) > ReadLen: #if full length and no primer found, exit, except if len is less than read length
                        continue
                #if full_length is passed, then only trim primers
                if not args.full_length:
                    #got here if primers were found they were trimmed
                    #now check seq length, pad if too short, trim if too long
                    if len(Seq) < args.min_len: #need this check here or primer dimers will get through
                        TooShort += 1
                        continue
                    if len(Seq) < args.trim_len and args.pad == 'on':
                        pad = args.trim_len - len(Seq)
                        Seq = Seq + pad*'N'
                        Qual = Qual + pad*'J'
                    else: #len(Seq) > args.trim_len:
                        Seq = Seq[:args.trim_len]
                        Qual = Qual[:args.trim_len]
                #got here, reads are primers trimmed and trim/padded, check length
                if len(Seq) < args.min_len:
                    TooShort += 1
                    continue
                ValidSeqs += 1     
                #now fix header
                header = indexReads.get(readID)
                Title = 'R_'+str(ValidSeqs)+';barcodelabel='+header[0]+';bcseq='+header[1]+';bcdiffs='+str(header[2])+';'
                #now write to file
                out.write("@%s\n%s\n+\n%s\n" % (Title, Seq, Qual))
            counts.write("%i,%i,%i,%i,%i,%i\n" % (Total, NoBC, NoPrimer, RevPrimerFound, TooShort, ValidSeqs))


args.out = re.sub(r'\W+', '', args.out)

log_name = args.out+'.amptk-demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)

amptklib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
amptklib.SystemInfo()
#get version of amptk
usearch = args.usearch
amptklib.versionDependencyChecks(usearch)

#get number of CPUs to use
if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

#parse a mapping file or a barcode fasta file, primers, etc get setup
#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
barcode_file = args.out + ".barcodes_used.fa"
if os.path.isfile(barcode_file):
    os.remove(barcode_file)

#check if mapping file passed, use this if present, otherwise use command line arguments
if args.mapping_file:
    if not os.path.isfile(args.mapping_file):
        amptklib.error("Mapping file is not valid: %s" % args.mapping_file)
        sys.exit(1)
    mapdata = amptklib.parseMappingFile(args.mapping_file, barcode_file)
    #forward primer in first item in tuple, reverse in second
    FwdPrimer = mapdata[0]
    RevPrimer = mapdata[1]
    genericmapfile = args.mapping_file
else:
    if args.barcode_fasta:
        if not os.path.isfile(args.barcode_fasta):
            amptklib.error("Mapping file or barcode_fasta is required")
            sys.exit(1)
        else:
            shutil.copyfile(args.barcode_fasta, barcode_file)
    
    #parse primers here so doesn't conflict with mapping primers
    #look up primer db otherwise default to entry
    if args.F_primer in amptklib.primer_db:
        FwdPrimer = amptklib.primer_db.get(args.F_primer)
    else:
        FwdPrimer = args.F_primer
    if args.R_primer in amptklib.primer_db:
        RevPrimer = amptklib.primer_db.get(args.R_primer)
    else:
        RevPrimer = args.R_primer

#setup 
if args.mapping_file:
    mapdict = amptklib.mapping2dict(args.mapping_file)
else:
    mapdict = False

amptklib.log.info("Loading %i samples from mapping file" % len(mapdict))

#process the index file, lookup in mapping file sample name, return dictionary
#will return dictionary:  readID : (SampleID, BC, mismatches) and list of reads to be discarded
amptklib.log.info("Mapping barcodes to sample IDs")
indexReads, discard = amptklib.barcodes2dict(args.index, mapdict, args.barcode_mismatch)

#estimate read length
if amptklib.check_valid_file(args.fastq):
    #if read length explicity passed use it otherwise measure it
    if args.read_length:
        ReadLen = args.read_length
    else:
        ReadLen = amptklib.GuessRL(args.fastq)

#create tmpdir
tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)

#Count FASTQ records
amptklib.log.info("Loading FASTQ Records")
orig_total = amptklib.countfastq(args.fastq)
size = amptklib.checkfastqsize(args.fastq)
readablesize = amptklib.convertSize(size)
amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#now we can merge the reads
mergedReads = args.out+'.merged.fastq'
amptklib.log.info("Merging PE reads using VSEARCH and filtering for Phix")
amptklib.MergeReads(args.fastq, args.reverse, tmpdir, mergedReads, ReadLen, args.min_len, args.usearch, args.rescue_forward, 'vsearch')

if not args.full_length:
    if args.pad == 'off':
        amptklib.log.info("Stripping primers and trim to %s bp" % (args.trim_len))
    else:
        amptklib.log.info("Stripping primers and trim/pad to %s bp" % (args.trim_len))
else:
    amptklib.log.info("Stripping primers and keeping only full length sequences")
amptklib.log.info("splitting the job over %i cpus, but this may still take awhile" % (cpus))

#now process the reads, have single file, so split like ion and run over multiple cores
#split fastq file
amptklib.split_fastq(os.path.join(tmpdir, mergedReads), orig_total, tmpdir, cpus*2)    


#start here to process the reads, first reverse complement the reverse primer
RevPrimer = revcomp_lib.RevComp(RevPrimer)
amptklib.log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

#now get file list from tmp folder
file_list = []
for file in os.listdir(tmpdir):
    if file.endswith(".fq"):
        file = os.path.join(tmpdir, file)
        file_list.append(file)

#finally process reads over number of cpus
amptklib.runMultiProgress(processRead, file_list, cpus)

print "-------------------------------------------------------"
#Now concatenate all of the demuxed files together
amptklib.log.info("Concatenating Demuxed Files")

Demux = args.out + '.demux.fq'
with open(Demux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(tmpdir,'*.demux.fq')):
        if filename == Demux:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)
#parse the stats
finalstats = [0,0,0,0,0,0]
for file in os.listdir(tmpdir):
    if file.endswith('.stats'):
        with open(os.path.join(tmpdir, file), 'rU') as statsfile:
            line = statsfile.readline()
            line = line.replace('\n', '')
            newstats = line.split(',')
            newstats = [int(i) for i in newstats]
            for x, num in enumerate(newstats):
                finalstats[x] += num
            
#clean up tmp folder
shutil.rmtree(tmpdir)

#output stats of the run
amptklib.log.info('{0:,}'.format(finalstats[0])+' total reads')
amptklib.log.info('{0:,}'.format(finalstats[1])+' discarded no index match')
amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1]-finalstats[2])+' Fwd Primer found, {0:,}'.format(finalstats[3])+ ' Rev Primer found')
amptklib.log.info('{0:,}'.format(finalstats[4])+' discarded too short (< %i bp)' % args.min_len)
amptklib.log.info('{0:,}'.format(finalstats[5])+' valid output reads')

#now loop through data and find barcoded samples, counting each.....
BarcodeCount = {}
with open(Demux, 'rU') as input:
    header = itertools.islice(input, 0, None, 4)
    for line in header:
        ID = line.split("=",1)[-1].split(";")[0]
        if ID not in BarcodeCount:
            BarcodeCount[ID] = 1
        else:
            BarcodeCount[ID] += 1

#now let's count the barcodes found and count the number of times they are found.
barcode_counts = "%30s:  %s" % ('Sample', 'Count')
for k,v in natsorted(BarcodeCount.items(), key=lambda (k,v): v, reverse=True):
    barcode_counts += "\n%30s:  %s" % (k, str(BarcodeCount[k]))
amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

#get file size
filesize = os.path.getsize(Demux)
readablesize = amptklib.convertSize(filesize)
amptklib.log.info("Output file:  %s (%s)" % (Demux, readablesize))
amptklib.log.info("Mapping file: %s" % args.mapping_file)
print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: amptk cluster -i %s -o out\n" % (Demux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (Demux)

sys.exit(1)


if args.barcode_fasta == 'pgm_barcodes.fa':
    amptklib.log.error("You did not specify a --barcode_fasta, it is required this type of data")
    os._exit(1)
if args.reverse:
    FNULL = open(os.devnull, 'w')
    #test for usearch
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        amptklib.log.warning("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    amptklib.log.info("USEARCH version: %s" % usearch_test)
    
    #next run USEARCH8 mergepe
    SeqIn = args.out + '.merged.fq'
    mergelog = args.out + '.mergedPE.log'
    notmerged = args.out + '.notmergedfwd.fq'
    amptklib.log.info("Merging PE Illumina reads with USEARCH")
    amptklib.log.debug("%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastq_truncqual 5 -fastq_maxdiffs 8 -minhsp 12" % (usearch, args.fastq, args.reverse, SeqIn))
    with open(mergelog, 'w') as logfile:
        subprocess.call([usearch, '-fastq_mergepairs', args.fastq, '-reverse', args.reverse, '-fastqout', SeqIn, '-fastq_truncqual', '5','-minhsp', '12','-fastq_maxdiffs', '8', '-fastqout_notmerged_fwd', notmerged], stdout = logfile, stderr = logfile)
    #recover forward reads that were not merged, could be longer sequence and still be ok
    with open(SeqIn, 'a') as output:
        with open(notmerged, 'rU') as input:
            for line in input:
                output.write(line)
else:
    amptklib.log.info("Running AMPtk on forward Illumina reads")
    SeqIn = args.fastq 

RevPrimer = revcomp_lib.RevComp(RevPrimer)
amptklib.log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
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

#setup barcode dictionary
Barcodes = fasta.ReadSeqsDict(barcode_file)

#setup for looking for reverse barcode
if args.reverse_barcode:
    RevBarcodes = {}
    rev_barcode_file = args.out + '.revbarcodes_used.fa'
    if os.path.isfile(rev_barcode_file):
        os.remove(rev_barcode_file)
    if not os.path.isfile(args.reverse_barcode):
        amptklib.log.info("Reverse barcode is not a valid file, exiting")
        sys.exit(1) 
    shutil.copyfile(args.reverse_barcode, rev_barcode_file)
    #parse and put into dictionary
    with open(rev_barcode_file, 'w') as output:
        with open(args.reverse_barcode, 'rU') as input:
            for rec in SeqIO.parse(input, 'fasta'):
                RevSeq = str(rec.seq.reverse_complement())
                if not RevSeq in RevBarcodes:
                    RevBarcodes[RevSeq] = rec.id
                else:
                    amptklib.log.error("Duplicate reverse barcodes detected, exiting")
                    sys.exit(1)
    
#get number of CPUs to use
if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

#get other values
MAX_PRIMER_MISMATCHES = int(args.primer_mismatch)
LabelPrefix = args.prefix
MinLen = int(args.min_len)
TrimLen = int(args.trim_len)
PL = len(FwdPrimer)
RL = len(RevPrimer)
OutCount = 0

#split the input FASTQ file into chunks to process
with open(SeqIn, 'rU') as input:
    SeqCount = amptklib.countfastq(SeqIn)
    amptklib.log.info('{0:,}'.format(SeqCount) + ' records loaded')
    SeqRecords = SeqIO.parse(SeqIn, 'fastq')
    chunks = SeqCount / cpus + 1
    amptklib.log.info("splitting job over %i cpus, this may still take awhile" % cpus)
    #divide into chunks, store in tmp file
    pid = os.getpid()
    folder = 'amptk_tmp_' + str(pid)
    if not os.path.exists(folder):
        os.makedirs(folder)
    for i, batch in enumerate(amptklib.batch_iterator(SeqRecords, chunks)) :
        filename = "chunk_%i.fq" % (i+1)
        tmpout = os.path.join(folder, filename)
        handle = open(tmpout, "w")
        count = SeqIO.write(batch, handle, "fastq")
        handle.close()

#now get file list from tmp folder
file_list = []
for file in os.listdir(folder):
    if file.endswith(".fq"):
        file = os.path.join(folder, file)
        file_list.append(file)

p = multiprocessing.Pool(cpus)
for f in file_list:
    p.apply_async(worker, [f])
p.close()
p.join()

print "-------------------------------------------------------"
#Now concatenate all of the demuxed files together
amptklib.log.info("Concatenating Demuxed Files")

tmpDemux = args.out + '.tmp.demux.fq'
with open(tmpDemux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(folder,'*.demux.fq')):
        if filename == tmpDemux:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)

#clean up tmp folder
shutil.rmtree(folder)

#last thing is to re-number of reads as it is possible they could have same name from multitprocessor split
catDemux = args.out + '.demux.fq'
amptklib.fastqreindex(tmpDemux, catDemux)
os.remove(tmpDemux)
        
amptklib.log.info("Counting FASTQ Records")
total = amptklib.countfastq(catDemux)
amptklib.log.info('{0:,}'.format(total) + ' reads processed')

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
amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

#get file size
filesize = os.path.getsize(catDemux)
readablesize = amptklib.convertSize(filesize)
amptklib.log.info("Output file:  %s (%s)" % (catDemux, readablesize))

print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: amptk cluster -i %s -o out\n" % (catDemux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (catDemux)
