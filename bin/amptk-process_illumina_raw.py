#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, logging, subprocess, multiprocessing, glob, itertools, re
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.revcomp_lib as revcomp_lib
import lib.amptklib as amptklib
import edlib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-process_illumina_raw.py', 
	usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--forward', dest='fastq', required=True, help='Illumina FASTQ R1 reads')
parser.add_argument('-r', '--reverse', required=True, help='Illumina FASTQ R2 reads')
parser.add_argument('-i', '--index', nargs='+', required=True, help='Illumina FASTQ index reads')
parser.add_argument('-m', '--mapping_file', help='QIIME-like mapping file')
parser.add_argument('--read_length', type=int, help='Read length, i.e. 2 x 300 bp = 300')
parser.add_argument('-o','--out', dest="out", default='illumina_out', help='Base name for output')
parser.add_argument('--fwd_primer', dest="F_primer", help='Forward Primer')
parser.add_argument('--rev_primer', dest="R_primer", help='Reverse Primer')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--barcode_mismatch', default=0, type=int, help='Number of mis-matches in barcode')
parser.add_argument('--barcode_fasta', help='FASTA file containing Barcodes (Names & Sequences)')
parser.add_argument('--require_primer', dest="primer", default='off', choices=['on', 'off'], help='Require Fwd primer to be present')
parser.add_argument('--rescue_forward', default='on', choices=['on', 'off'], help='Rescue Not-merged forward reads')
parser.add_argument('--barcode_rev_comp', action='store_true', help='Reverse complement barcode sequences')
parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default=300, type=int, help='Trim length for reads')
parser.add_argument('-p','--pad', default='off', choices=['on', 'off'], help='Pad with Ns to a set length')
parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
parser.add_argument('--cleanup', action='store_true', help='remove intermediate files')
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
                #set postions to zero
                ForCutPos, RevCutPos = (None,)*2
                #first thing is look for forward primer, if found trim it off
                ForCutPos = amptklib.findFwdPrimer(FwdPrimer, seq, args.primer_mismatch, amptklib.degenNuc)
                #if require primer is on make finding primer in amplicon required if amplicon is larger than read length
                #if less than read length, can't enforce primer because could have been trimmed via staggered trim in fastq_mergepairs
                if args.primer == 'on' and len(seq) > ReadLen:
                    if not ForCutPos:
                        NoPrimer += 1
                        continue
                    Seq = seq[ForCutPos:]
                    Qual = qual[ForCutPos:]
                else:
                    if ForCutPos:
                        Seq = seq[ForCutPos:]
                        Qual = qual[ForCutPos:]
                    else:
                        NoPrimer += 1
                        Seq = seq
                        Qual = qual
                #now look for reverse primer
                RevCutPos = amptklib.findRevPrimer(RevPrimer, Seq, args.primer_mismatch, amptklib.degenNuc)
                if RevCutPos:
                    RevPrimerFound += 1
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
                #now write to file
                out.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))
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

#create tmpdir
tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)

#parse a mapping file or a barcode fasta file, primers, etc get setup
#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
barcode_file = args.out + ".barcodes_used.fa"
if os.path.isfile(barcode_file):
    os.remove(barcode_file)

#check if mapping file passed, use this if present, otherwise use command line arguments
if args.mapping_file:
    if not os.path.isfile(args.mapping_file):
        amptklib.log.error("Mapping file is not valid: %s" % args.mapping_file)
        sys.exit(1)
    mapdata = amptklib.parseMappingFile(args.mapping_file, barcode_file)
    #forward primer in first item in tuple, reverse in second
    FwdPrimer = mapdata[0]
    RevPrimer = mapdata[1]
    genericmapfile = args.mapping_file
		
if not FwdPrimer or not RevPrimer:
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

#if still no primers set, then exit
if not FwdPrimer or not RevPrimer:
    amptklib.log.error("Please provide primer sequences via --fwd_primer and --rev_primer")
    sys.exit(1)

#setup 
if args.mapping_file:
    mapdict = amptklib.mapping2dict(args.mapping_file)
else:
	if not args.barcode_fasta:
		amptklib.log.error("A -m,--mapping_file or --barcode_fasta are required")
		sys.exit(1)
	else:
    	mapdict = {}
    	with open(args.barcode_fasta, 'rU') as infile:
    		for rec in SeqIO.parse(infile, 'fasta'):
    			Seq = str(rec.seq)
    			if not Seq in mapdict:
    				mapdict[Seq] = rec.id

#if barcodes_rev_comp passed then reverse complement the keys in mapdict
if args.barcode_rev_comp:
	amptklib.log.info("Reverse complementing barcode sequences")
	backupDict = mapdict
	mapdict = {}
	for k,v in backupDict.items():
		RCkey = amptklib.RevComp(k)
		if not RCkey in mapdict:
			mapdict[RCkey] = v

amptklib.log.info("Loading %i samples from mapping file, checking FASTQ input" % len(mapdict))

#rename reads according to indexes
if not amptklib.PEandIndexCheck(args.fastq, args.reverse, args.index[0]): #check they are all same length
    amptklib.log.error("FASTQ input malformed, read numbers do not match")
    sys.exit(1)
amptklib.log.info("Mapping indexes to reads and renaming PE reads")
cleanR1 = os.path.join(tmpdir, 'renamedR1.fastq')
cleanR2 = os.path.join(tmpdir, 'renamedR2.fastq')
amptklib.DemuxIllumina(args.fastq, args.reverse, args.index[0], mapdict, args.barcode_mismatch, cleanR1, cleanR2)

amptklib.log.info("Loading FASTQ Records")
#estimate read length
if amptklib.check_valid_file(cleanR1):
    #if read length explicity passed use it otherwise measure it
    if args.read_length:
        ReadLen = args.read_length
    else:
        ReadLen = amptklib.GuessRL(cleanR1)
        amptklib.log.info('Estimation of read length is %i bp' % ReadLen)

#Count FASTQ records
orig_total = amptklib.countfastq(cleanR1)
size = amptklib.checkfastqsize(cleanR1)
readablesize = amptklib.convertSize(size)
amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#now we can merge the reads
mergedReads = args.out+'.merged.fastq'
amptklib.log.info("Merging PE reads using VSEARCH and filtering for phiX")
amptklib.MergeReads(os.path.abspath(cleanR1), os.path.abspath(cleanR2), tmpdir, mergedReads, ReadLen, args.min_len, args.usearch, args.rescue_forward, 'vsearch', '', 1)

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
mergedTotal = amptklib.countfastq(os.path.join(tmpdir, mergedReads))
amptklib.split_fastq(os.path.join(tmpdir, mergedReads), mergedTotal, tmpdir, cpus*2)    


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

tmpDemux = os.path.join(tmpdir, args.out + '.demux.fq')
with open(tmpDemux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(tmpdir,'*.demux.fq')):
        if filename == tmpDemux:
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

#finally reindex output
#last thing is to re-number of reads as it is possible they could have same name from multitprocessor split
Demux = args.out + '.demux.fq'
amptklib.fastqreindex(tmpDemux, Demux)
os.remove(tmpDemux)

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

#create mapping file if one doesn't exist
genericmapfile = args.out + '.mapping_file.txt'
amptklib.CreateGenericMappingFile(barcode_file, FwdPrimer, revcomp_lib.RevComp(RevPrimer), '', genericmapfile, BarcodeCount)

#compress the output to save space
FinalDemux = Demux+'.gz'
amptklib.Fzip(Demux, FinalDemux, cpus)
amptklib.removefile(Demux)

if args.cleanup:
	amptklib.SafeRemove(tmpdir)

#get file size
filesize = os.path.getsize(FinalDemux)
readablesize = amptklib.convertSize(filesize)
amptklib.log.info("Output file:  %s (%s)" % (FinalDemux, readablesize))
amptklib.log.info("Mapping file: %s" % args.mapping_file)
print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: amptk cluster -i %s -o out\n" % (FinalDemux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (FinalDemux)
