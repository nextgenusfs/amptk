#!/usr/bin/env python

#This script runs UPARSE_REF and then generates some stats
#written by Jon Palmer palmer.jona at gmail dot com

import sys
import os
import argparse
import subprocess
import csv

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

parser=argparse.ArgumentParser(prog='ficus-uparse_ref_summary.py',
    description='''Script runs uparse_ref and then generates counts from the output per each member of mock community. ''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('fastq', help='fastq')
parser.add_argument('db', help='DB of mock community')
parser.add_argument('-o','--out', dest='out', default="out", help='Basename for output files')
parser.add_argument('--cluster_otus', dest='otus', action='store_true', help='Run OTU clustering before UPARSE_REF')
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

if args.otus:
    #run clustering prior to Uparse-ref
    #now run usearch8 full length dereplication
    print bcolors.BLUE + "Running Dereplication" + bcolors.ENDC
    print "------------------------------------------------"
    derep_out = args.out + '.derep.fq'
    os.system('%s %s %s %s %s %s' % (usearch, '-derep_fulllength', args.fastq, '-sizeout', '-fastqout', derep_out))

    #now run usearch 8 sort by size
    print "------------------------------------------------"
    print bcolors.BLUE + "Running SortBySize" + bcolors.ENDC
    print "------------------------------------------------"
    sort_out = args.out + '.sort.fa'
    os.system('%s %s %s %s %s'% (usearch, '-sortbysize', derep_out, '-minsize 2 -fastaout', sort_out))
    #now run clustering algorithm
    print "------------------------------------------------"
    print bcolors.BLUE + "Running UPARSE Clustering at 97 percent" + bcolors.ENDC
    print "------------------------------------------------"
    otu_out = args.out + '.otus.fa'
    os.system('%s %s %s %s %s %s %s %s' % (usearch, '-cluster_otus', sort_out, '-sizein -sizeout -relabel OTU_', '-otu_radius_pct', '3', '-otus', otu_out))
    uparse_input = otu_out
    
else:
    print bcolors.BLUE + "Running Dereplication" + bcolors.ENDC
    print "------------------------------------------------"
    derep_out = args.out + '.derep.fq'
    os.system('%s %s %s %s %s %s' % (usearch, '-derep_fulllength', args.fastq, '-sizeout', '-fastqout', derep_out))
    uparse_input = derep_out

#now run usearch8 uparse_ref
print "------------------------------------------------"
print bcolors.BLUE + "Running UPARSE_ref" + bcolors.ENDC
print "------------------------------------------------"
uparse_out = args.out + '.uparse.txt'
os.system('%s %s %s %s %s %s %s' % (usearch, '-uparse_ref', uparse_input, '-strand plus -db', args.db, '-uparseout', uparse_out))
print "------------------------------------------------"
print bcolors.BLUE + "Analyzing Results" + bcolors.ENDC
#now try to parse results and get summary for each member of mock community
delete_list = ["barcodelabel=", "size=", ""]
input = open(uparse_out)
temp_out = open('temp.uparse.out', "w+")
for line in input:
    for word in delete_list:
        line = line.replace(word, "")
    line = line.replace(";", "\t")
    line = line.replace("\t\t", "\t")
    temp_out.write(line)
input.close()
temp_out.close()  

#get the mock community header names
db_file = open(args.db)
headerList = []
for line in db_file:
    if line[0] == ">":
        headerList.append(line[1:].strip())
db_file.close()
headerList = sorted(headerList)

#now read in csv and do some counting of reads and write to file
out_filename = args.out + ".uparse.stats.txt"
if os.path.isfile(out_filename):
    os.remove(out_filename)
out_file = open(out_filename, "a")
first_col = "Mock_ID"
head = "Perfect\tGood\tNoisy\tOther\tTotal"
out_file.write ("%s\t%s\n" % (first_col, head))
total_count = 0
for header in headerList:
    perfect_count = 0
    good_count = 0
    other_count = 0
    noisy_count = 0
    uparse_in = csv.reader(open('temp.uparse.out'), delimiter='\t')
    for row in uparse_in:
        if row[-1] == header and row[-4] == "perfect":
            perfect_count = perfect_count + int(row[-5])
        if row[-1] == header and row[-4] == "good":
            good_count = good_count + int(row[-5])
        if row[-1] == header and row[-4] == "other":
            other_count = other_count + int(row[-5])
        if row[-1] == header and row[-4] == "noisy":
            noisy_count = noisy_count + int(row[-5])
        row_total = perfect_count + good_count + other_count + noisy_count
    total_count = total_count + row_total
    out_file.write ("%s\t%i\t%i\t%i\t%i\t%i\n" % (header, perfect_count, good_count, noisy_count, other_count, row_total))
chimera_count = 0
nohit_count = 0
uparse_in2 = csv.reader(open('temp.uparse.out'), delimiter='\t')
for row2 in uparse_in2:
    if row2[-4] == "chimera":
        chimera_count = chimera_count + int(row2[-5])
    if row2[-1] == "*":
        nohit_count = nohit_count + int(row2[-5])
raw_total = chimera_count + total_count + nohit_count
print "Total Reads: " + str(raw_total)
print "NonChimeras: " + str(total_count)
print "Chimeras: " + str(chimera_count)
print "No Hits: " + str(nohit_count)
out_file.write ("Total Reads\t%i\n" % (raw_total))
out_file.write ("NonChimeras\t%i\n" % (total_count))
out_file.write ("Chimeric Reads\t%i\n" % (chimera_count))
out_file.write ("No Hit\t%i" % (nohit_count))
os.remove('temp.uparse.out')
out_file.close()
print bcolors.GREEN + "Done!" + bcolors.ENDC
print "------------------------------------------------"