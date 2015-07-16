#!/usr/bin/env python

#This script runs USEARCH OTU clustering
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import subprocess
import inspect

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='ficus-OTU_cluster.py',
    description='''Script runs UPARSE OTU clustering. 
    Requires USEARCH and uc2otutab.py by Robert Edgar: http://drive5.com''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('fastq', help='FASTQ file from ficus-process.py')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
parser.add_argument('-m','--minsize', default='2', help='Min size to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--uchime_ref', default='False', choices=['ITS1','ITS2','Full'], help='Run UCHIME REF')
parser.add_argument('--map_unfiltered', action='store_true', help='map original reads back to OTUs')

args=parser.parse_args()

#open log file for usearch8 stderr redirect
log_name = args.out + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)
log_file = open(log_name, 'ab')

usearch = args.usearch
try:
    subprocess.call([usearch, '--version'], stdout = log_file, stderr = log_file)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)
        
#now run usearch8 fastq filtering step
filter_out = args.out + '.EE' + args.maxee + '.filter.fq'
print "\nCMD: %s -fastq_filter %s -fastq_maxee %s -fastqout %s\n" % (usearch, args.fastq, args.maxee, filter_out)
subprocess.call([usearch, '-fastq_filter', args.fastq, '-fastq_maxee', args.maxee, '-fastqout', filter_out], stdout = log_file, stderr = log_file)

#now run usearch8 full length dereplication
derep_out = args.out + '.EE' + args.maxee + '.derep.fa'
print "CMD: %s -derep_fulllength %s -sizeout -fastaout %s\n" % (usearch, filter_out, derep_out)
subprocess.call([usearch, '-derep_fulllength', filter_out, '-sizeout', '-fastaout', derep_out], stdout = log_file, stderr = log_file)

#now run usearch 8 sort by size
sort_out = args.out + '.EE' + args.maxee + '.sort.fa'
print "CMD: %s -sortbysize %s -minsize %s -fastaout %s\n" % (usearch, derep_out, args.minsize, sort_out)
subprocess.call([usearch, '-sortbysize', derep_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = log_file, stderr = log_file)

#now run clustering algorithm
radius = str(100 - int(args.pct_otu))
otu_out = args.out + '.EE' + args.maxee + '.otus.fa'
print "CMD: %s -cluster_otus %s -sizein -sizeout -relabel OTU_ -otu_radius_pct %s -otus %s\n" % (usearch, sort_out, radius, otu_out)
subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-sizeout', '-relabel', 'OTU_', '-otu_radius_pct', radius, '-otus', otu_out], stdout = log_file, stderr = log_file)

#optional UCHIME Ref 
if args.uchime_ref == "False":
    uchime_out = otu_out
else:
    uchime_out = args.out + '.EE' + args.maxee + '.uchime.otus.fa'
    #You might need to update these in the future, but leaving data and version in name so it is obvious where they came from
    if args.uchime_ref == "ITS1":
        uchime_db = script_path + "/lib/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta"
    if args.uchime_ref == "ITS2":
        uchime_db = script_path + "/lib/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta"
    if args.uchime_ref == "Full":
        uchime_db = script_path + "/lib/uchime_sh_refs_dynamic_original_985_11.03.2015.fasta"
    print "CMD: %s -uchime_ref %s -strand plus -db %s -nonchimeras %s\n" % (usearch, otu_out, uchime_db, uchime_out)
    subprocess.call([usearch, '-uchime_ref', otu_out, '-strand', 'plus', '-db', uchime_db, '-nonchimeras', uchime_out], stdout = log_file, stderr = log_file)
    
#now map reads back to OTUs
uc_out = args.out + '.EE' + args.maxee + '.mapping.uc'
if args.map_unfiltered:
    reads = args.fastq
else:
    reads = filter_out
print "CMD: %s -usearch_global %s -strand plus -id 0.97 -db %s -uc %s\n" % (usearch, reads, uchime_out, uc_out)
subprocess.call([usearch, '-usearch_global', reads, '-strand', 'plus', '-id', '0.97', '-db', uchime_out, '-uc', uc_out], stdout = log_file, stderr = log_file)
#Build OTU table
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"
print "CMD: python %s %s > %s" % (uc2tab, uc_out, otu_table)
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input FASTQ:           %s" % (args.fastq))
print ("Filtered FASTQ:        %s" % (filter_out))
print ("Dereplicated FASTA:    %s" % (derep_out))
print ("Sorted FASTA:          %s" % (sort_out))
print ("Clustered OTUs:        %s" % (otu_out))
if args.uchime_ref != "False":
    print ("Chimera Filtered OTUs: %s" % (uchime_out))
print ("UCLUST Mapping file:   %s" % (uc_out))
print ("OTU Table:             %s" % (otu_table))
print ("USEARCH LogFile:       %s" % (log_name))
print "------------------------------------------------"
