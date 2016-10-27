#!/usr/bin/env python

import sys, os, argparse, logging, shutil, subprocess, inspect
from Bio.SeqIO.QualityIO import FastqGeneralIterator
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-dada2.py',
    description='''Script takes output from ufits pre-processing and runs DADA2''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', required=True, help='Input Demuxed containing FASTQ')
parser.add_argument('-o','--out', default='dada2', help='Output Basename')
parser.add_argument('-l','--length', type=int, required=True, help='Length to truncate reads')
parser.add_argument('-m','--maxee', default='1.0', help='MaxEE quality filtering')
parser.add_argument('-p','--platform', default='ion', choices=['ion', 'illumina', '454'], help='Sequencing platform')
parser.add_argument('--cleanup', action='store_true', help='Remove Intermediate Files')
args=parser.parse_args()

dada2script = os.path.join(parentdir, 'bin', 'dada2_pipeline_nofilt.R')

def folder2list(input, ending):
    names = []
    if not os.path.isdir(input):
        return False
    else:
        for x in os.listdir(input):
            if x.endswith(ending):
                x = os.path.join(input, x)
                names.append(x)
    return names

def maxEE(input, maxee, output):
    subprocess.call(['vsearch', '--fastq_filter', input, '--fastq_maxee', str(maxee), '--fastqout', output, '--fastq_qmax', '45'], stdout=FNULL, stderr = FNULL)

def splitDemux(input, outputdir, length):
    for title, seq, qual in FastqGeneralIterator(open(input)):
        sample = title.split('barcodelabel=')[1]
        sample = sample.replace(';', '')
        Seq = seq.rstrip('N')
        Qual = qual[:len(Seq)]
        if len(Seq) >= int(length):
            with open(os.path.join(outputdir, sample+'.fastq'), 'ab') as output:
                output.write("@%s\n%s\n+\n%s\n" % (title, Seq[:int(length):], Qual[:int(length)]))

def getAvgLength(input):
    AvgLength = []
    for title, seq, qual in FastqGeneralIterator(open(input)):
        Seq = seq.rstrip('N')
        AvgLength.append(len(Seq))
    Average = sum(AvgLength) / float(len(AvgLength))
    Min = min(AvgLength)
    Max = max(AvgLength)
    return (Average, Min, Max)

#remove logfile if exists
log_name = args.out + '.ufits-dada2.log'
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
version = ufitslib.get_version()
ufitslib.log.info("%s" % version)

#check dependencies
programs = ['Rscript']
ufitslib.CheckDependencies(programs)

#check if vsearch version > 1.9.1 is installed
vsearch_check = ufitslib.which('vsearch')
if vsearch_check:
    vsearch = ufitslib.checkvsearch()
    vsearch_version = ufitslib.get_vsearch_version()
    if vsearch:
        ufitslib.log.info("VSEARCH v%s" % vsearch_version)
    else:
        ufitslib.log.info("VSEARCH v%s detected, need version at least v1.9.1, using Python for filtering")
else:
    vsearch = False
    ufitslib.log.info("VSEARCH not installed, using Python for filtering")    

#Count FASTQ records
ufitslib.log.info("Loading FASTQ Records")
orig_total = ufitslib.countfastq(args.fastq)
size = ufitslib.checkfastqsize(args.fastq)
readablesize = ufitslib.convertSize(size)
ufitslib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#quality filter
ufitslib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
derep = args.out+'.derep.fq'
if vsearch:
    maxEE(args.fastq, float(args.maxee), derep)
else:
    with open(derep, 'w') as output:
        SeqIO.write(ufitslib.MaxEEFilter(args.fastq, args.maxee), output, 'fastq')

total = ufitslib.countfastq(derep)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#Get Average length without any N's
averageLen = getAvgLength(derep)
ufitslib.log.info("DADA2 compatible read lengths, average: %i bp, minimum: %i bp, maximum: %i" % (averageLen[0], averageLen[1], averageLen[2]))
if averageLen[0] < int(args.length):
    TruncLen = int(averageLen[0] - 1)
    ufitslib.log.error('Warning: Average length of reads %i bp, is less than specified truncation length %s bp' % (averageLen[0], args.length))
    ufitslib.log.error('Resetting truncation length to %i bp' % TruncLen)
else:
    TruncLen = int(args.length)

#now split into individual files
ufitslib.log.info("Splitting FASTQ file by Sample and truncating to %i bp" % TruncLen)
filtfolder = args.out+'_filtered'
if not os.path.isdir(filtfolder):
    os.makedirs(filtfolder)
splitDemux(derep, filtfolder, TruncLen)

#now run DADA2 on filtered folder
ufitslib.log.info("Running DADA2 pipeline")
dada2log = args.out+'.dada2.R.log'
dada2out = args.out+'.dada2.csv'
with open(dada2log, 'w') as logfile:
    subprocess.call(['Rscript', '--vanilla', dada2script, filtfolder, dada2out, args.platform], stdout = logfile, stderr = logfile)

#check for results
if not os.path.isfile(dada2out):
    ufitslib.log.error("DADA2 run failed, please check %s logfile" % dada2log)
    sys.exit(1)
    
#now process the output, pull out fasta, rename, etc
ufitslib.log.info("Reformating DADA2 output")
fastaout = args.out+'.otus.fa'
otutable = args.out+'.otu_table.txt'
counter = 0
with open(fastaout, 'w') as writefasta:
    with open(otutable, 'w') as writetab:
        with open(dada2out, 'rU') as input:
            for line in input:
                line = line.replace('\n', '')
                line = line.replace('"', '')
                cols = line.split(',')
                if counter == 0:
                    header = "OTUId\t"+ '\t'.join(cols[1:]) + '\n'
                    writetab.write(header)
                else:
                    Seq = cols[0]
                    ID = 'OTU_'+str(counter)
                    newline = ID+'\t'+'\t'.join(cols[1:]) + '\n'
                    writetab.write(newline)
                    writefasta.write(">%s\n%s\n" % (ID, Seq))
                counter += 1

#get number of bimeras from logfile
with open(dada2log, 'rU') as bimeracheck:
    for line in bimeracheck:
        if line.startswith('Identified '):
            bimeraline = line.split(' ')
            bimeras = int(bimeraline[1])
            totalSeqs = int(bimeraline[5])
validSeqs = totalSeqs - bimeras
ufitslib.log.info('{0:,}'.format(totalSeqs) + ' total sequences')
ufitslib.log.info('{0:,}'.format(bimeras) + ' chimeras removed')
ufitslib.log.info('{0:,}'.format(validSeqs) + ' valid sequences')

if args.cleanup:
    shutil.rmtree(filtfolder)
    os.remove(dada2out)

#Print location of files to STDOUT
print "-------------------------------------------------------"
print "DADA2 Script has Finished Successfully"
print "-------------------------------------------------------"
if not args.cleanup:
    print "Tmp Folder of files: %s" % filtfolder
print "Clustered OTUs: %s" % fastaout
print "OTU Table: %s" % otutable
print "-------------------------------------------------------"

otu_print = fastaout.split('/')[-1]
tab_print = otutable.split('/')[-1]
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)

        