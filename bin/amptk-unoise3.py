#!/usr/bin/env python

#This script runs USEARCH UNOISE3 denoising
#written by Jon Palmer nextgenusfs@gmail.com

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import sys
import os
import argparse
import subprocess
import inspect
import logging
import shutil
import multiprocessing
from Bio import SeqIO
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.amptklib as amptklib


#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(script_path)

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-unoise2.py', usage="%(prog)s [options] -i file.demux.fq\n%(prog)s -h for help menu",
    description='''Script runs UNOISE2 algorithm.
    Requires USEARCH9 by Robert C. Edgar: http://drive5.com/usearch''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-m','--minsize', default='8', help='Min size to keep for denoising')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch10', help='USEARCH10 EXE')
parser.add_argument('-p','--pct_otu', default='97', help="Biological OTU Clustering Percent")
parser.add_argument('--uchime_ref', help='Run UCHIME2 REF [ITS,16S,LSU,COI,custom]')
parser.add_argument('--map_filtered', action='store_true', help='map quality filtered reads back to OTUs')
parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
args=parser.parse_args()

def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize

#remove logfile if exists
log_name = args.out + '.amptk-unoise3.log'
if os.path.isfile(log_name):
    os.remove(log_name)

amptklib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print("-------------------------------------------------------")

#initialize script, log system info and usearch version
amptklib.SystemInfo()
#Do a version check
usearch = args.usearch
amptklib.checkusearch10(usearch)

#make tmp folder
tmp = args.out + '_tmp'
if not os.path.exists(tmp):
    os.makedirs(tmp)

#Count FASTQ records
amptklib.log.info("Loading FASTQ Records")
#convert to FASTA for mapping
orig_fasta = os.path.join(tmp, args.out+'.orig.fa')
cmd = ['vsearch', '--fastq_filter', args.FASTQ, '--fastaout', orig_fasta, '--fastq_qmax', '55']
amptklib.runSubprocess(cmd, amptklib.log)
orig_total = amptklib.countfasta(orig_fasta)
size = amptklib.checkfastqsize(args.FASTQ)
readablesize = amptklib.convertSize(size)
amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#Expected Errors filtering step
filter_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.filter.fq')
filter_fasta = os.path.join(tmp, args.out + '.EE' + args.maxee + '.filter.fa')
amptklib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
cmd = ['vsearch', '--fastq_filter', args.FASTQ, '--fastq_maxee', str(args.maxee), '--fastqout', filter_out, '--fastaout', filter_fasta, '--fastq_qmax', '55']
amptklib.runSubprocess(cmd, amptklib.log)
total = amptklib.countfastq(filter_out)
amptklib.log.info('{0:,}'.format(total) + ' reads passed')

#now run full length dereplication
derep_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.derep.fa')
amptklib.log.info("De-replication (remove duplicate reads)")
cmd = ['vsearch', '--derep_fulllength', filter_out, '--relabel', 'Read_', '--sizeout', '--output', derep_out]
amptklib.runSubprocess(cmd, amptklib.log)
total = amptklib.countfasta(derep_out)
amptklib.log.info('{0:,}'.format(total) + ' reads passed')

#now run de-noiser UNOISE3
amptklib.log.info("Denoising reads with UNOISE3")
unoise_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.unoise.fa')
cmd = [usearch, '-unoise3', derep_out, '-zotus', unoise_out, '-minsize', args.minsize]
amptklib.runSubprocess(cmd, amptklib.log)
total = amptklib.countfasta(unoise_out)
amptklib.log.info('{0:,}'.format(total) + ' denoised sequences')

#strip N's
#amptklib.log.info("Cleaning up padding from OTUs")
otu_clean = os.path.join(tmp, args.out + '.EE' + args.maxee + '.clean.fa')
amptklib.fasta_strip_padding(unoise_out, otu_clean)

#run optional uchime_ref
if not args.uchime_ref:
    uchime_out = otu_clean
else:
    uchime_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.uchime.otus.fa')
    #R. Edgar now says using largest DB is better for UCHIME, so use the one distributed with taxonomy
    if args.uchime_ref in ['ITS', '16S', 'LSU', 'COI']: #test if it is one that is setup, otherwise default to full path
        uchime_db = os.path.join(parentdir, 'DB', args.uchime_ref+'.extracted.fa')
        if not os.path.isfile(uchime_db):
            amptklib.log.error("Database not properly configured, run `amptk install` to setup DB, skipping chimera filtering")
            uchime_out = otu_clean
    else:
        uchime_db = os.path.abspath(args.uchime_ref)
    #now run chimera filtering if all checks out
    if not os.path.isfile(uchime_out):
        amptklib.log.info("Chimera Filtering (VSEARCH)")
        cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', otu_clean, '--db', uchime_db, '--nonchimeras', uchime_out]
        amptklib.runSubprocess(cmd, amptklib.log)
        total = amptklib.countfasta(uchime_out)
        amptklib.log.info('{0:,}'.format(total) + ' OTUs passed')

#inferred sequences
iSeqs = args.out+'.iSeqs.fa'
amptklib.fastarename(uchime_out, 'iSeq', iSeqs)

#build OTU table with iSeqs
uc_iSeq_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.mapping.uc')
iSeq_otu_table = args.out + '.otu_table.txt'
#setup reads to map
if args.map_filtered:
    reads = filter_fasta
else:
    reads = orig_fasta
amptklib.log.info("Mapping Reads to iSeqs and Building OTU table")
cmd = ['vsearch', '--usearch_global', reads, '--strand', 'plus', '--id', '0.97', '--db', iSeqs, '--uc', uc_iSeq_out, '--otutabout', iSeq_otu_table]
amptklib.runSubprocess(cmd, amptklib.log)

#count reads mapped
total = amptklib.line_count2(uc_iSeq_out)
amptklib.log.info('{0:,}'.format(total) + ' reads mapped to iSeqs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

#now cluster to biological OTUs with UCLUST
radius = float(args.pct_otu) / 100.
amptklib.log.info("Clustering denoised sequences into biological OTUs at %s%%" % args.pct_otu)
uclust_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.uclust.fa')
cmd = ['vsearch', '--cluster_smallmem', iSeqs, '--centroids', uclust_out, '--id', str(radius), '--strand', 'plus', '--relabel', 'OTU', '--qmask', 'none', '--usersort']
amptklib.runSubprocess(cmd, amptklib.log)
total = amptklib.countfasta(uclust_out)
amptklib.log.info('{0:,}'.format(total) + ' OTUs generated')

#determine where denoised sequences clustered
ClusterComp = args.out+'.iSeqs2clusters.txt'
iSeqmap = args.out+'.unoise_map.uc'
cmd = [usearch, '-usearch_global', iSeqs, '-db', uclust_out, '-id', str(radius), '-uc', iSeqmap, '-strand', 'plus']
amptklib.runSubprocess(cmd, amptklib.log)
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
    for k,v in natsorted(list(iSeqMapped.items())):
        clusters.write('%s\t%s\n' % (k, ', '.join(v)))

#now map reads back to OTUs and build OTU table
uc_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.cluster.mapping.uc')
otu_table = os.path.join(tmp, args.out + '.EE' + args.maxee + '.cluster.otu_table.txt')
#setup reads to map
if args.map_filtered:
    reads = filter_fasta
else:
    reads = orig_fasta
amptklib.log.info("Mapping Reads to OTUs and Building OTU table")
cmd = ['vsearch', '--usearch_global', reads, '--strand', 'plus', '--id', '0.97', '--db', uclust_out, '--uc', uc_out, '--otutabout', otu_table]
amptklib.runSubprocess(cmd, amptklib.log)

#count reads mapped
total = amptklib.line_count2(uc_out)
amptklib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

#Move files around, delete tmp if argument passed.
currentdir = os.getcwd()
final_otu = os.path.join(currentdir, args.out + '.cluster.otus.fa')
shutil.copyfile(uclust_out, final_otu)
final_otu_table = os.path.join(currentdir, args.out + '.cluster.otu_table.txt')
shutil.copyfile(otu_table, final_otu_table)
if not args.debug:
    shutil.rmtree(tmp)

#Print location of files to STDOUT
print("-------------------------------------------------------")
print("UNOISE2 Script has Finished Successfully")
print("-------------------------------------------------------")
if not not args.debug:
    print("Tmp Folder of files: %s" % tmp)
print("inferred Seqs: %s" % os.path.abspath(iSeqs))
print("iSeq OTU Table: %s" % os.path.abspath(iSeq_otu_table))
print("Clustered OTUs: %s" % final_otu)
print("OTU Table: %s" % final_otu_table)
print("iSeqs 2 OTUs: %s" % os.path.abspath(ClusterComp))
print("-------------------------------------------------------")

otu_print = final_otu.split('/')[-1]
tab_print = final_otu_table.split('/')[-1]
if 'win32' in sys.platform:
    print("\nExample of next cmd: amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print))
else:
    print(colr.WARN + "\nExample of next cmd:" + colr.END + " amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print))
