#!/usr/bin/env python

#This script runs USEARCH UNOISE2 denoising
#written by Jon Palmer nextgenusfs@gmail.com

import sys, os, argparse, subprocess, inspect, csv, re, logging, shutil, multiprocessing
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib
from natsort import natsorted

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(script_path)

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-unoise2.py', usage="%(prog)s [options] -i file.demux.fq\n%(prog)s -h for help menu",
    description='''Script runs UNOISE2 algorithm.
    Requires USEARCH9 by Robert C. Edgar: http://drive5.com/usearch''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-m','--minampout', default='4', help='Min size to keep for denoising')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
parser.add_argument('-p','--pct_otu', default='97', help="Biological OTU Clustering Percent")
parser.add_argument('--uchime_ref', help='Run UCHIME2 REF [ITS,16S,LSU,COI,custom]')
parser.add_argument('--map_filtered', action='store_true', help='map quality filtered reads back to OTUs')
parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
args=parser.parse_args()

def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize

#remove logfile if exists
log_name = args.out + '.ufits-unoise2.log'
if os.path.isfile(log_name):
    os.remove(log_name)

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

#make tmp folder
tmp = args.out + '_tmp'
if not os.path.exists(tmp):
    os.makedirs(tmp)

#Count FASTQ records
ufitslib.log.info("Loading FASTQ Records")
orig_total = ufitslib.countfastq(args.FASTQ)
size = checkfastqsize(args.FASTQ)
readablesize = ufitslib.convertSize(size)
ufitslib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#Expected Errors filtering step and convert to fasta
filter_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.filter.fq')
filter_fasta = os.path.join(tmp, args.out + '.EE' + args.maxee + '.filter.fa')
orig_fasta = os.path.join(tmp, args.out+'.orig.fa')
ufitslib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
cmd = ['vsearch', '--fastq_filter', args.FASTQ, '--fastq_maxee', str(args.maxee), '--fastqout', filter_out, '--fastaout', filter_fasta, '--fastq_qmax', '55']
ufitslib.runSubprocess(cmd, ufitslib.log)
cmd = ['vsearch', '--fastq_filter', args.FASTQ, '--fastaout', orig_fasta, '--fastq_qmax', '55']
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.countfastq(filter_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#now run full length dereplication
derep_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.derep.fa')
ufitslib.log.info("De-replication (remove duplicate reads)")
cmd = ['vsearch', '--derep_fulllength', filter_out, '--relabel', 'Read_', '--sizeout', '--output', derep_out]
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.countfasta(derep_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#now run de-noiser UNOISE2
ufitslib.log.info("Denoising reads with UNOISE2")
unoise_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.unoise.fa')
cmd = [usearch, '-unoise2', derep_out, '-fastaout', unoise_out, '--minampsize', args.minampout]
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.countfasta(unoise_out)
ufitslib.log.info('{0:,}'.format(total) + ' denoised sequences')

#now cluster to biological OTUs with UCLUST
radius = float(args.pct_otu) / 100.
ufitslib.log.info("Clustering denoised sequences into OTUs at %s%%" % args.pct_otu)
uclust_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.uclust.fa')
cmd = [usearch, '-cluster_smallmem', unoise_out, '-id', str(radius), '-centroids', uclust_out, '-relabel', 'OTU']
ufitslib.runSubprocess(cmd, ufitslib.log)
total = ufitslib.countfasta(uclust_out)
ufitslib.log.info('{0:,}'.format(total) + ' OTUs generated')

#determine where denoised sequences clustered
ClusterComp = args.out+'.denoised2clusters.txt'
iSeqmap = args.out+'.unoise_map.uc'
cmd = [usearch, '-usearch_global', unoise_out, '-db', uclust_out, '-id', str(radius), '-uc', iSeqmap, '-strand', 'plus']
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

#strip N's
ufitslib.log.info("Cleaning up padding from OTUs")
otu_clean = os.path.join(tmp, args.out + '.EE' + args.maxee + '.clean.fa')
ufitslib.fasta_strip_padding(uclust_out, otu_clean)


#run optional uchime_ref
if not args.uchime_ref:
    uchime_out = otu_clean
else:
    uchime_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.uchime.otus.fa')
    #R. Edgar now says using largest DB is better for UCHIME, so use the one distributed with taxonomy
    if args.uchime_ref in ['ITS', '16S', 'LSU', 'COI']: #test if it is one that is setup, otherwise default to full path
        uchime_db = os.path.join(parentdir, 'DB', args.uchime_ref+'.extracted.fa')
        if not os.path.isfile(uchime_db):
            ufitslib.log.error("Database not properly configured, run `ufits install` to setup DB, skipping chimera filtering")
            uchime_out = otu_clean
    else:
        uchime_db = os.path.abspath(args.uchime_ref)
    #now run chimera filtering if all checks out
    if not os.path.isfile(uchime_out):
        ufitslib.log.info("Chimera Filtering (VSEARCH)")
        cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', otu_clean, '--db', uchime_db, '--nonchimeras', uchime_out]
        ufitslib.runSubprocess(cmd, ufitslib.log)
        total = ufitslib.countfasta(uchime_out)
        ufitslib.log.info('{0:,}'.format(total) + ' OTUs passed')

#now map reads back to OTUs and build OTU table
uc_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.mapping.uc')
otu_table = os.path.join(tmp, args.out + '.EE' + args.maxee + '.otu_table.txt')
#setup reads to map
if args.map_filtered:
    reads = filter_fasta
else:
    reads = orig_fasta
ufitslib.log.info("Mapping Reads to OTUs and Building OTU table")
cmd = ['vsearch', '--usearch_global', reads, '--strand', 'plus', '--id', '0.97', '--db', uchime_out, '--uc', uc_out, '--otutabout', otu_table]
ufitslib.runSubprocess(cmd, ufitslib.log)

#count reads mapped
total = ufitslib.line_count(uc_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

#Move files around, delete tmp if argument passed.
currentdir = os.getcwd()
final_otu = os.path.join(currentdir, args.out + '.cluster.otus.fa')
shutil.copyfile(uchime_out, final_otu)
final_otu_table = os.path.join(currentdir, args.out + '.otu_table.txt')
shutil.copyfile(otu_table, final_otu_table)
if not args.debug:
    shutil.rmtree(tmp)

#Print location of files to STDOUT
print "-------------------------------------------------------"
print "UNOISE2 Script has Finished Successfully"
print "-------------------------------------------------------"
if not not args.debug:
    print "Tmp Folder of files: %s" % tmp
print "Clustered OTUs: %s" % final_otu
print "OTU Table: %s" % final_otu_table
print "-------------------------------------------------------"

otu_print = final_otu.split('/')[-1]
tab_print = final_otu_table.split('/')[-1]
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
