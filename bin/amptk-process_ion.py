#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, logging, subprocess, multiprocessing, glob, itertools, re, gzip
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

parser=argparse.ArgumentParser(prog='amptk-process_ion.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
    description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq','--sff', '--fasta', '--bam', dest='fastq', required=True, help='BAM/FASTQ/SFF/FASTA file')
parser.add_argument('-q','--qual', help='QUAL file (if -i is FASTA)')
parser.add_argument('-o','--out', dest="out", default='ion', help='Base name for output')
parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer')
parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer')
parser.add_argument('-m','--mapping_file', help='Mapping file: QIIME format can have extra meta data columns')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
parser.add_argument('--barcode_fasta', default='pgm_barcodes.fa', help='FASTA file containing Barcodes (Names & Sequences)')
parser.add_argument('--reverse_barcode', help='FASTA file containing 3 prime Barocdes')
parser.add_argument('-b','--list_barcodes', dest="barcodes", default='all', help='Enter Barcodes used separated by commas')
parser.add_argument('--min_len', default=50, type=int, help='Minimum read length to keep')
parser.add_argument('-l','--trim_len', default=250, type=int, help='Trim length for reads')
parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
parser.add_argument('--mult_samples', dest="multi", default='False', help='Combine multiple samples (i.e. FACE1)')
parser.add_argument('--illumina', action='store_true', help='Input data is single file Illumina')
parser.add_argument('--ion', action='store_true', help='Input data is Ion Torrent')
parser.add_argument('--454', action='store_true', help='Input data is 454')
parser.add_argument('--reverse', help='Illumina reverse reads')
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH EXE')
args=parser.parse_args()

def FindBarcode(Seq, BarcodeDict):
    for BarcodeLabel in BarcodeDict.keys():
        Barcode = BarcodeDict[BarcodeLabel]
        if Seq.startswith(Barcode):
            return Barcode, BarcodeLabel
    return "", ""

def processRead(input):
    base = os.path.basename(input).split('.')[0]
    PL = len(FwdPrimer)
    RL = len(RevPrimer)
    with open(os.path.join(tmpdir, base+'.demux.fq'), 'w') as out:
        counter = 0
        for title, seq, qual in FastqGeneralIterator(open(input)):
            bcmismatch = 0
            #look for barcode
            Barcode, BarcodeLabel = FindBarcode(seq, Barcodes)
            if Barcode == "":
                continue
            BarcodeLength = len(Barcode)
            #now search for forward primer
            BestPosFor, BestDiffsFor = primer.BestMatch2(seq, FwdPrimer, args.primer_mismatch)
            if BestPosFor > 0 and BestPosFor <= BarcodeLength+2: #if found will be > 0, and should be found after barcode
                ForTrim = BestPosFor + PL
                #now search for reverse primer
                BestPosRev, BestDiffsRev = primer.BestMatch2(seq, RevPrimer, args.primer_mismatch)
                if BestPosRev > 0:  #reverse primer was found    
                    #location to trim sequences
                    RevTrim = BestPosRev                
                    #determine reverse barcode
                    if args.reverse_barcode:
                        RevBCdiffs = 0
                        BCcut = BestPosRev + RL
                        CutSeq = seq[BCcut:]
                        RevBarcode, RevBarcodeLabel = FindBarcode(CutSeq, RevBarcodes)
                        if RevBarcode == "":
                            continue
                        BarcodeLabel = BarcodeLabel+'_'+RevBarcodeLabel                       
                    #now trim record remove forward and reverse reads
                    Seq = seq[ForTrim:RevTrim]
                    Qual = qual[ForTrim:RevTrim]
                    #since found reverse primer, now also need to pad/trim
                    if not args.full_length:
                        if len(Seq) < args.trim_len:
                            pad = args.trim_len - len(Seq)
                            Seq = Seq + pad*'N'
                            Qual = Qual +pad*'J'
                        elif len(Seq) > args.trim_len:
                            Seq = Seq[:args.trim_len]
                            Qual = Qual[:args.trim_len]
                else:
                    #trim record, did not find reverse primer
                    if args.full_length: #if full length then move to next record
                        continue
                    #trim away forward primer
                    Seq = seq[ForTrim:]
                    Qual = seq[ForTrim:]
                    #check length and trim, throw away if too short as it was bad read
                    if len(Seq) < args.trim_len:
                        continue
                    Seq = Seq[:args.trim_len]
                    Qual = Qual[:args.trim_len]
                #check minimum length
                if len(Seq) >= int(args.min_len):
                    counter += 1
                    #rename header
                    Name = 'R_'+str(counter)+';barcodelabel='+BarcodeLabel+';'
                    out.write("@%s\n%s\n+\n%s\n" % (Name, Seq, Qual))

args.out = re.sub(r'\W+', '', args.out)

log_name = args.out + '.amptk-demux.log'
if os.path.isfile(log_name):
    os.remove(log_name)
FNULL = open(os.devnull, 'w')
amptklib.setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
amptklib.SystemInfo()
#Do a version check
usearch = args.usearch
amptklib.versionDependencyChecks(usearch)

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
    if args.barcode_fasta == 'pgm_barcodes.fa':
        #get script path and barcode file name
        pgm_barcodes = os.path.join(parentdir, 'DB', args.barcode_fasta)
        if args.barcodes == "all":
            if args.multi == 'False':
                shutil.copyfile(pgm_barcodes, barcode_file)
            else:
                with open(barcode_file, 'w') as barcodeout:
                    with open(pgm_barcodes, 'rU') as input:
                        for rec in SeqIO.parse(input, 'fasta'):
                            outname = args.multi+'.'+rec.id
                            barcodeout.write(">%s\n%s\n" % (outname, rec.seq))
        else:
            bc_list = args.barcodes.split(",")
            inputSeqFile = open(pgm_barcodes, "rU")
            SeqRecords = SeqIO.to_dict(SeqIO.parse(inputSeqFile, "fasta"))
            for rec in bc_list:
                name = "BC." + rec
                seq = SeqRecords[name].seq
                if args.multi != 'False':
                    outname = args.multi+'.'+name
                else:
                    outname = name
                outputSeqFile = open(barcode_file, "a")
                outputSeqFile.write(">%s\n%s\n" % (outname, seq))
            outputSeqFile.close()
            inputSeqFile.close()
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

    #because we use an 'A' linker between barcode and primer sequence, add an A if ion is chemistry
    if args.ion:
        FwdPrimer = 'A' + FwdPrimer
        Adapter = 'CCATCTCATCCCTGCGTGTCTCCGACTCAG'
    else:
        Adapter = ''

#check if input is compressed
gzip_list = []
if args.fastq.endswith('.gz'):
    gzip_list.append(os.path.abspath(args.fastq))
if args.reverse:
    if args.reverse.endswith('.gz'):
        gzip_list.append(os.path.abspath(args.reverse))
if gzip_list:
    amptklib.log.info("Gzipped input files detected, uncompressing")
    for file in gzip_list:
        amptklib.log.debug("Uncompressing %s" % file)
        OutName = os.path.splitext(file)[0]
        InFile = gzip.open(file, 'rU')
        ReadFile = InFile.read()
        OutFile = open(OutName, 'w')
        OutFile.write(ReadFile)
        OutFile.close()
        InFile.close()
        os.remove(file) #remove .gz file  
    args.fastq = args.fastq.replace('.gz', '')
    if args.reverse:
        args.reverse = args.reverse.replace('.gz', '')
        
#if SFF file passed, convert to FASTQ with biopython
if args.fastq.endswith(".sff"):
    if args.barcode_fasta == 'pgm_barcodes.fa':
        if not args.mapping_file:
            amptklib.log.error("You did not specify a --barcode_fasta or --mapping_file, one is required for 454 data")
            sys.exit(1)
    amptklib.log.info("SFF input detected, converting to FASTQ")
    SeqIn = args.out + '.sff.extract.fastq'
    SeqIO.convert(args.fastq, "sff-trim", SeqIn, "fastq")
elif args.fastq.endswith(".fas") or args.fastq.endswith(".fasta") or args.fastq.endswith(".fa"):
    if not args.qual:
        amptklib.log.error("FASTA input detected, however no QUAL file was given.  You must have FASTA + QUAL files")
        sys.exit(1)
    else:
        if args.barcode_fasta == 'pgm_barcodes.fa':
            if not args.mapping_file:
                amptklib.log.error("You did not specify a --barcode_fasta or --mapping_file, one is required for 454 data")
                sys.exit(1)
        SeqIn = args.out + '.fastq'
        amptklib.log.info("FASTA + QUAL detected, converting to FASTQ")
        amptklib.faqual2fastq(args.fastq, args.qual, SeqIn)
elif args.fastq.endswith('.bam'):
    amptklib.CheckDependencies(['bedtools'])
    SeqIn = args.out+'.fastq'
    amptklib.log.info("Converting Ion Torrent BAM file to FASTQ using BedTools")
    cmd = ['bedtools', 'bamtofastq', '-i', args.fastq, '-fq', SeqIn]
    amptklib.runSubprocess(cmd, amptklib.log)
else:        
    SeqIn = args.fastq

#check if illumina argument is passed, if so then run merge PE
if args.illumina:
    if args.barcode_fasta == 'pgm_barcodes.fa':
        if not args.mapping_file:
            amptklib.log.error("You did not specify a --barcode_fasta or --mapping_file, one is required for Illumina2 data")
            sys.exit(1)
    if args.reverse:
        #next run USEARCH9 mergePE
        #get read length
        RL = amptklib.GuessRL(args.fastq)
        #merge reads
        amptklib.log.info("Merging Illumina reads")
        SeqIn = args.out + '.merged.fq'
        amptklib.MergeReads(args.fastq, args.reverse, '.', SeqIn, RL, args.min_len, usearch, 'on')
    else:
        amptklib.log.info("Running AMPtk on forward Illumina reads")
        SeqIn = args.fastq 

#start here to process the reads, first reverse complement the reverse primer
RevPrimer = revcomp_lib.RevComp(RevPrimer)
amptklib.log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

#then setup barcode dictionary
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
                if not rec.id in RevBarcodes:
                    RevBarcodes[rec.id] = RevSeq
                    output.write('>%s\n%s\n' % (rec.id, RevSeq))
                else:
                    amptklib.log.error("Duplicate reverse barcodes detected, exiting")
                    sys.exit(1)
#get number of CPUs to use
if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

#Count FASTQ records
amptklib.log.info("Loading FASTQ Records")
orig_total = amptklib.countfastq(SeqIn)
size = amptklib.checkfastqsize(SeqIn)
readablesize = amptklib.convertSize(size)
amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#create tmpdir and split input into n cpus
tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
#split the input FASTQ file into chunks to process
with open(SeqIn, 'rU') as input:
    SeqRecords = SeqIO.parse(SeqIn, 'fastq')
    chunks = orig_total / (2*cpus)+1
    #divide into chunks, store in tmp file
    for i, batch in enumerate(amptklib.batch_iterator(SeqRecords, chunks)) :
        filename = "chunk_%i.fq" % (i+1)
        tmpout = os.path.join(tmpdir, filename)
        handle = open(tmpout, "w")
        count = SeqIO.write(batch, handle, "fastq")
        handle.close()

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

tmpDemux = args.out + '.tmp.demux.fq'
with open(tmpDemux, 'wb') as outfile:
    for filename in glob.glob(os.path.join(tmpdir,'*.demux.fq')):
        if filename == tmpDemux:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)

#clean up tmp folder
shutil.rmtree(tmpdir)

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
        ID = line.split("=",1)[-1].split(";")[0]
        if ID not in BarcodeCount:
            BarcodeCount[ID] = 1
        else:
            BarcodeCount[ID] += 1

#now let's count the barcodes found and count the number of times they are found.
barcode_counts = "%30s:  %s" % ('Sample', 'Count')
barcodes_found = []
for k,v in natsorted(BarcodeCount.items(), key=lambda (k,v): v, reverse=True):
    barcode_counts += "\n%30s:  %s" % (k, str(BarcodeCount[k]))
    if k not in barcodes_found:
        barcodes_found.append(k)
amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

if not args.mapping_file:
    #create a generic mappingfile for downstream processes
    genericmapfile = args.out + '.mapping_file.txt'
    amptklib.CreateGenericMappingFile(barcode_file, FwdPrimer, revcomp_lib.RevComp(RevPrimer), Adapter, genericmapfile, barcodes_found)

#get file size
filesize = os.path.getsize(catDemux)
readablesize = amptklib.convertSize(filesize)
amptklib.log.info("Output file:  %s (%s)" % (catDemux, readablesize))
amptklib.log.info("Mapping file: %s" % genericmapfile)

print "-------------------------------------------------------"
if 'win32' in sys.platform:
    print "\nExample of next cmd: amptk cluster -i %s -o out\n" % (catDemux)
else:
    print col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (catDemux)
