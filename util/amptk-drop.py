#!/usr/bin/env python

#This script renames fasta
import sys, argparse, os, inspect
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.amptklib as amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser = argparse.ArgumentParser(prog='amptk-drop.py',
    description = '''Script that drops OTUs and then creates OTU table''',
    epilog = """Written by Jon Palmer (2016) nextgenusfs@gmail.com""", formatter_class=MyFormatter)
parser.add_argument('-i', '--input', required=True, help='OTUs in FASTA format')
parser.add_argument('-r', '--reads', required=True, help='Demuxed reads FASTQ format')
parser.add_argument('-o','--out', default='amptk-drop', help='Base output name')
parser.add_argument('-l','--list', nargs='+', help='Input list of (BC) names to remove')
parser.add_argument('-f','--file', help='File containing list of names to remove')
args=parser.parse_args()

#remove logfile if exists
log_name = args.out + '.amptk-drop.log'
if os.path.isfile(log_name):
    os.remove(log_name)

amptklib.setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
amptklib.SystemInfo()

#check the list or file parameters, one of them must have something
if not args.list:
    if not args.file:
        amptklib.log.error("Error, you must specifiy a list of OTU names or a file containing names")
        sys.exit(1)
if not args.file:
    if not args.list:
        amptklib.log.error("Error, you must specifiy a list of OTU names or a file containing names")
        sys.exit(1)
if args.list and args.file:
    amptklib.log.error("Error, you must specifiy either list of OTU names or a file containing OTU names, not both")
    sys.exit(1)
if args.file:   
    count = amptklib.line_count(args.file)
    #load in list of names to remove
    with open(args.file, 'rU') as input:
        lines = [line.rstrip('\n') for line in input]
if args.list:
    count = len(args.list)
    lines = args.list
#make sure it is a set, faster lookup
dropList = set(lines)

#load data
total = amptklib.countfasta(args.input)
amptklib.log.info("Loading %i OTUs" % total)

#load in the fasta file, change if in dictionary and output to stdout
amptklib.log.info("Dropping %i OTUs" % count)
newOTUs = args.out + '.cleaned.otus.fa'
with open(newOTUs, 'w') as otus:
    with open(args.input, 'rU') as fasta:
        for rec in SeqIO.parse(fasta, 'fasta'):
            if not rec.id in dropList:
                SeqIO.write(rec, otus, 'fasta')

#now make new OTU table
amptklib.log.info("Mapping Reads to OTUs and Building OTU table")
newTable = args.out + '.cleaned.otu_table.txt'
tmpReads = args.out + '.reads.tmp'
uc_out = args.out + '.mapping.uc'
cmd = ['vsearch', '--fastq_filter', args.reads, '--fastaout', tmpReads, '--fastq_qmax', '55']
amptklib.runSubprocess(cmd, amptklib.log)
cmd = ['vsearch', '--usearch_global', tmpReads, '--strand', 'plus', '--id', '0.97', '--db', newOTUs, '--uc', uc_out, '--otutabout', newTable]
amptklib.runSubprocess(cmd, amptklib.log)

#count OTUs
otu_count = amptklib.countfasta(newOTUs)
amptklib.log.info('{0:,}'.format(otu_count) + ' OTUs remaining')

#count reads mapped
total = amptklib.line_count(uc_out)
orig_total = amptklib.countfasta(tmpReads)
amptklib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

#Print location of files to STDOUT
print "-------------------------------------------------------"
print "Clustered OTUs: %s" % newOTUs
print "OTU Table: %s" % newTable
print "-------------------------------------------------------"

#cleanup
amptklib.removefile(tmpReads)
amptklib.removefile(uc_out)

otu_print = newOTUs.split('/')[-1]
tab_print = newTable.split('/')[-1]
if 'win32' in sys.platform:
    print "\nExample of next cmd: amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)


