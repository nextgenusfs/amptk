#!/usr/bin/env python

#This script runs USEARCH OTU clustering
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import subprocess
import inspect

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class bcolors:
    GREEN = '\033[92m'
    BLUE = '\033[36m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)           

parser=argparse.ArgumentParser(prog='ficus-OTU_cluster.py',
    description='''Script runs UPARSE OTU clustering. 
    Requires USEARCH and uc2otutab.py by Robert Edgar: http://drive5.com''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('fastq', help='FASTQ file from fastq_strip_relabel.py')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
parser.add_argument('--uchime_ref', default='False', choices=['ITS1','ITS2','Full'], help='Run UCHIME REF')
parser.add_argument('--keep_singletons', action='store_true', help='Keep singletons before clustering')
parser.add_argument('--map_filtered_reads', action='store_true', help='map quality trimmed reads back to OTUs')
args=parser.parse_args()

usearch = "usearch8"
try:
    print "------------------------------------------------"
    print bcolors.BLUE + "Looking for USEARCH8 in your PATH:" + bcolors.ENDC
    subprocess.call([usearch, '--version'])
    print "------------------------------------------------"
except OSError:
    print + bcolors.FAIL + "%s not found in your PATH" % usearch + bcolors.ENDC
    usearch = raw_input("Enter full path to USEARCH or type exit to quit: ")
    
if usearch == "exit":
    os._exit(1)
    
#now run usearch8 fastq filtering step
print bcolors.BLUE + "Running FASTQ Filtering" + bcolors.ENDC
print "------------------------------------------------"
filter_out = args.out + '.EE' + args.maxee + '.filter.fq'
os.system('%s %s %s %s %s %s %s' % (usearch, '-fastq_filter', args.fastq, '-fastq_maxee', args.maxee, '-fastqout', filter_out))

#now run usearch8 full length dereplication
print "------------------------------------------------"
print bcolors.BLUE + "Running Dereplication" + bcolors.ENDC
print "------------------------------------------------"
derep_out = args.out + '.EE' + args.maxee + '.derep.fa'
os.system('%s %s %s %s %s %s' % (usearch, '-derep_fulllength', filter_out, '-sizeout', '-fastaout', derep_out))

#now run usearch 8 sort by size
print "------------------------------------------------"
print bcolors.BLUE + "Running SortBySize" + bcolors.ENDC
print "------------------------------------------------"
sort_out = args.out + '.EE' + args.maxee + '.sort.fa'
if args.keep_singletons:
    singletons = "1"
else:
    singletons = "2"
os.system('%s %s %s %s %s %s %s' % (usearch, '-sortbysize', derep_out, '-minsize', singletons, '-fastaout', sort_out))

#now run clustering algorithm
radius = str(100 - int(args.pct_otu))
print "------------------------------------------------"
print bcolors.BLUE + "Running UPARSE Clustering at %s percent" % (args.pct_otu) + bcolors.ENDC
print "------------------------------------------------"
otu_out = args.out + '.EE' + args.maxee + '.otus.fa'
os.system('%s %s %s %s %s %s %s %s' % (usearch, '-cluster_otus', sort_out, '-sizein -sizeout -relabel OTU_', '-otu_radius_pct', radius, '-otus', otu_out))

#optional UCHIME Ref 
if args.uchime_ref == "False":
    uchime_out = otu_out
else:
    print "------------------------------------------------"
    print bcolors.BLUE + "Running UCHIME-Ref" + bcolors.ENDC
    print "------------------------------------------------"
    uchime_out = args.out + '.EE' + args.maxee + '.uchime.fa'
    #You might need to update these in the future, but leaving data and version in name so it is obvious where they came from
    if args.uchime_ref == "ITS1":
        uchime_db = script_path + "/lib/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta"
    if args.uchime_ref == "ITS2":
        uchime_db = script_path + "/lib/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta"
    if args.uchime_ref == "Full":
        uchime_db = script_path + "/lib/uchime_sh_refs_dynamic_original_985_11.03.2015.fasta"
    os.system('%s %s %s %s %s %s %s' % (usearch, '-uchime_ref', otu_out, '-strand plus -db', uchime_db, '-nonchimeras', uchime_out))
    
#now map reads back to OTUs
print "------------------------------------------------"
print bcolors.BLUE + "Mapping Reads to OTUs with usearch_global" + bcolors.ENDC
print "------------------------------------------------"
uc_out = args.out + '.EE' + args.maxee + '.mapping.uc'
if args.map_filtered_reads:
    reads = filter_out
else:
    reads = args.fastq
os.system('%s %s %s %s %s %s %s' % (usearch, '-usearch_global', reads, '-strand plus -id 0.97 -db', uchime_out, '-uc', uc_out))
#Build OTU table
print "------------------------------------------------"
print bcolors.BLUE + "Converting to OTU table" + bcolors.ENDC
print "------------------------------------------------"
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))
print "------------------------------------------------"
print bcolors.GREEN + "OTU Clustering Script has Finished Successfully" + bcolors.ENDC
print "------------------------------------------------"