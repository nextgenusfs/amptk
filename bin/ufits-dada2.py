#!/usr/bin/env python

import sys, os, argparse, logging, shutil, subprocess, inspect
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib
import numpy as np
from natsort import natsorted

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
parser.add_argument('-e','--maxee', default='1.0', help='MaxEE quality filtering')
parser.add_argument('-p','--pct_otu', default='97', help="Biological OTU Clustering Percent")
parser.add_argument('--platform', default='ion', choices=['ion', 'illumina', '454'], help='Sequencing platform')
parser.add_argument('--uchime_ref', help='Run UCHIME REF [ITS,16S,LSU,COI,custom]')
parser.add_argument('--pool', action='store_true', help='Pool all sequences together for DADA2')
parser.add_argument('--debug', action='store_true', help='Keep all intermediate files')
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

def splitDemux(input, outputdir, length):
    for title, seq, qual in FastqGeneralIterator(open(input)):
        sample = title.split('barcodelabel=')[1]
        sample = sample.replace(';', '')
        if len(seq) >= int(length):
            with open(os.path.join(outputdir, sample+'.fastq'), 'ab') as output:
                output.write("@%s\n%s\n+\n%s\n" % (title, seq[:int(length):], qual[:int(length)]))

def getAvgLength(input):
    AvgLength = []
    for title, seq, qual in FastqGeneralIterator(open(input)):
        AvgLength.append(len(seq))
    Average = sum(AvgLength) / float(len(AvgLength))
    Min = min(AvgLength)
    Max = max(AvgLength)
    a = np.array(AvgLength)
    nintyfive = np.percentile(a, 5)
    return (Average, Min, Max, int(nintyfive))

#remove logfile if exists
log_name = args.out + '.ufits-dada2.log'
if os.path.isfile(log_name):
    ufitslib.removefile(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"
#initialize script, log system info and usearch version
ufitslib.SystemInfo()
#Do a version check
usearch = args.usearch
ufitslib.versionDependencyChecks(usearch)

#check dependencies
programs = ['Rscript']
ufitslib.CheckDependencies(programs)

#Count FASTQ records and remove 3' N's as dada2 can't handle them
ufitslib.log.info("Loading FASTQ Records")
orig_total = ufitslib.countfastq(args.fastq)
size = ufitslib.checkfastqsize(args.fastq)
readablesize = ufitslib.convertSize(size)
ufitslib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')
no_ns = args.out+'.cleaned_input.fq'
ufitslib.fastq_strip_padding(args.fastq, no_ns)

#quality filter
ufitslib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
derep = args.out+'.qual-filtered.fq'
filtercmd = ['vsearch', '--fastq_filter', no_ns, '--fastq_maxee', str(args.maxee), '--fastqout', derep, '--fastq_qmax', '55', '--fastq_maxns', '0']
ufitslib.runSubprocess(filtercmd, ufitslib.log)
total = ufitslib.countfastq(derep)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#Get Average length without any N's
averageLen = getAvgLength(derep)
ufitslib.log.info("DADA2 compatible read lengths, avg: %i bp, min: %i bp, max: %i bp, top 95%%: %i bp" % (averageLen[0], averageLen[1], averageLen[2], averageLen[3]))
if averageLen[0] < int(args.length):
    TruncLen = int(averageLen[3])
    ufitslib.log.error('Warning: Average length of reads %i bp, is less than specified truncation length %s bp' % (averageLen[0], args.length))
    ufitslib.log.error('Resetting truncation length to %i bp (keep > 95%% of data) ' % TruncLen)
else:
    TruncLen = int(args.length)

#now split into individual files
ufitslib.log.info("Splitting FASTQ file by Sample and truncating to %i bp" % TruncLen)
filtfolder = args.out+'_filtered'
if os.path.isdir(filtfolder):
    shutil.rmtree(filtfolder)
os.makedirs(filtfolder)
splitDemux(derep, filtfolder, TruncLen)

#now run DADA2 on filtered folder
ufitslib.log.info("Running DADA2 pipeline")
dada2log = args.out+'.dada2.Rscript.log'
dada2out = args.out+'.dada2.csv'
#check pooling vs notpooled, default is not pooled.
if args.pool:
    POOL = 'TRUE'
else:
    POOL = 'FALSE'
CORES = str(ufitslib.getCPUS())
with open(dada2log, 'w') as logfile:
    subprocess.call(['Rscript', '--vanilla', dada2script, filtfolder, dada2out, args.platform, POOL, CORES], stdout = logfile, stderr = logfile)

#check for results
if not os.path.isfile(dada2out):
    ufitslib.log.error("DADA2 run failed, please check %s logfile" % dada2log)
    sys.exit(1)
    
#now process the output, pull out fasta, rename, etc
fastaout = args.out+'.otus.tmp'
counter = 1
with open(fastaout, 'w') as writefasta:
    with open(dada2out, 'rU') as input:
        next(input)
        for line in input:
            line = line.replace('\n', '')
            line = line.replace('"', '')
            cols = line.split(',')
            Seq = cols[0]
            ID = 'iSeq_'+str(counter)
            writefasta.write(">%s\n%s\n" % (ID, Seq))
            counter += 1

#get number of bimeras from logfile
with open(dada2log, 'rU') as bimeracheck:
    for line in bimeracheck:
        if line.startswith('Identified '):
            bimeraline = line.split(' ')
            bimeras = int(bimeraline[1])
            totalSeqs = int(bimeraline[5])
        if line.startswith('[1] "dada2'):
            dada2version = line.split(' ')[-1].replace('"\n', '').rstrip()
        if line.startswith('[1] "R '):
            Rversion = line.split(' ')[-1].replace('"\n', '').rstrip()
validSeqs = totalSeqs - bimeras
ufitslib.log.info("R v%s, DADA2 v%s" % (Rversion, dada2version))
ufitslib.log.info('{0:,}'.format(totalSeqs) + ' total inferred sequences (iSeqs)')
ufitslib.log.info('{0:,}'.format(bimeras) + ' denovo chimeras removed')
ufitslib.log.info('{0:,}'.format(validSeqs) + ' valid iSeqs')

#optional UCHIME Ref
uchime_out = args.out+'.nonchimeras.fa'
chimeraFreeTable = args.out+'.otu_table.txt'
iSeqs = args.out+'.iSeqs.fa'
if not args.uchime_ref:
    os.rename(fastaout, iSeqs)
else:
    #check if file is present, remove from previous run if it is.
    if os.path.isfile(uchime_out):
        ufitslib.removefile(uchime_out)
    #R. Edgar now says using largest DB is better for UCHIME, so use the one distributed with taxonomy
    if args.uchime_ref in ['ITS', '16S', 'LSU', 'COI']: #test if it is one that is setup, otherwise default to full path
        uchime_db = os.path.join(parentdir, 'DB', args.uchime_ref+'.extracted.fa')
        if not os.path.isfile(uchime_db):
            ufitslib.log.error("Database not properly configured, run `ufits install` to setup DB, skipping chimera filtering")
            uchime_out = fastaout
    else:
        if os.path.isfile(args.uchime_ref):
            uchime_db = os.path.abspath(args.uchime_ref)
        else:
            ufitslib.log.error("%s is not a valid file, skipping reference chimera filtering" % args.uchime_ref)
            uchime_out = fastaout
    #now run chimera filtering if all checks out
    if not os.path.isfile(uchime_out):
        ufitslib.log.info("Chimera Filtering (VSEARCH) using %s DB" % args.uchime_ref)
        cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', fastaout, '--db', uchime_db, '--nonchimeras', uchime_out]
        ufitslib.runSubprocess(cmd, ufitslib.log)
        total = ufitslib.countfasta(uchime_out)
        uchime_chimeras = validSeqs - total
        ufitslib.log.info('{0:,}'.format(total) + ' iSeqs passed, ' + '{0:,}'.format(uchime_chimeras) + ' ref chimeras removed')

    #now reformat OTUs and OTU table, dropping chimeric OTUs from table, sorting the output as well
    nonchimeras = ufitslib.fasta2list(uchime_out)
    inferredSeqs = SeqIO.index(uchime_out, 'fasta')
    with open(iSeqs, 'w') as iSeqout:
        for x in natsorted(nonchimeras):
            SeqIO.write(inferredSeqs[x], iSeqout, 'fasta')
    if not args.debug:
        #clean up chimeras fasta
        ufitslib.removefile(uchime_out)
        if os.path.isfile(fastaout):
            ufitslib.removefile(fastaout)


#setup output files
dadademux = args.out+'.dada2.map.uc'
bioSeqs = args.out+'.cluster.otus.fa'
bioTable = args.out+'.cluster.otu_table.txt'
demuxtmp = args.out+'.original.fa'
uctmp = args.out+'.map.uc'
ClusterComp = args.out+'.iSeqs2clusters.txt'

#map reads to DADA2 OTUs
ufitslib.log.info("Mapping reads to DADA2 iSeqs")
cmd = ['vsearch', '--fastq_filter', os.path.abspath(args.fastq),'--fastq_qmax', '55', '--fastq_maxns', '0', '--fastaout', demuxtmp]
ufitslib.runSubprocess(cmd, ufitslib.log)
cmd = ['vsearch', '--usearch_global', demuxtmp, '--db', iSeqs, '--id', '0.97', '--uc', dadademux, '--strand', 'plus', '--otutabout', chimeraFreeTable ]
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.line_count(dadademux)
ufitslib.log.info('{0:,}'.format(total) + ' reads mapped to iSeqs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

#cluster
ufitslib.log.info("Clustering iSeqs at %s%% to generate biological OTUs" % args.pct_otu)
radius = float(args.pct_otu) / 100.
cmd = ['vsearch', '--cluster_smallmem', iSeqs, '--centroids', bioSeqs, '--id', str(radius), '--strand', 'plus', '--relabel', 'OTU', '--qmask', 'none', '--usersort']
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.countfasta(bioSeqs)
ufitslib.log.info('{0:,}'.format(total) + ' OTUs generated')

#determine where iSeqs clustered
iSeqmap = args.out+'.iseq_map.uc'
cmd = ['vsearch', '--usearch_global', iSeqs, '--db', bioSeqs, '--id', str(radius), '--uc', iSeqmap, '--strand', 'plus']
ufitslib.runSubprocess(cmd, ufitslib.log)
iSeqMapped = {}
with open(iSeqmap, 'rU') as mapping:
    for line in mapping:
        line = line.replace('\n', '')
        cols = line.split('\t')
        OTU = cols[9]
        Hit = cols[8]
        if not OTU in iSeqMapped:
            iSeqMapped[OTU] = [Hit]
        else:
            iSeqMapped[OTU].append(Hit)
with open(ClusterComp, 'w') as clusters:
    clusters.write('OTU\tiSeqs\n')
    for k,v in natsorted(iSeqMapped.items()):
        clusters.write('%s\t%s\n' % (k, ', '.join(v)))
#create OTU table
ufitslib.log.info("Mapping reads to OTUs")
cmd = ['vsearch', '--usearch_global', demuxtmp, '--db', bioSeqs, '--id', '0.97', '--uc', uctmp, '--strand', 'plus', '--otutabout', bioTable]
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.line_count(uctmp)
ufitslib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

if not args.debug:
    ufitslib.removefile(no_ns)
    shutil.rmtree(filtfolder)
    ufitslib.removefile(dada2out)
    ufitslib.removefile(derep)
    ufitslib.removefile(demuxtmp)
    ufitslib.removefile(uctmp)
    ufitslib.removefile(iSeqmap)
    ufitslib.removefile(dadademux)

#Print location of files to STDOUT
print "-------------------------------------------------------"
print "DADA2 Script has Finished Successfully"
print "-------------------------------------------------------"
if args.debug:
    print "Tmp Folder of files: %s" % filtfolder
print "Inferred iSeqs: %s" % iSeqs
print "iSeq OTU Table: %s" % chimeraFreeTable
print "Clustered OTUs: %s" % bioSeqs
print "OTU Table: %s" % bioTable
print "iSeqs 2 OTUs: %s" % ClusterComp
print "-------------------------------------------------------"

otu_print = bioSeqs.split('/')[-1]
tab_print = bioTable.split('/')[-1]
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)

        