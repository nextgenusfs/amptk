#!/usr/bin/env python

#This script runs USEARCH OTU clustering
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import subprocess
import inspect
import csv
import re

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='ufits-OTU_cluster.py', usage="%(prog)s [options] -f file.demux.fq\n%(prog)s -h for help menu",
    description='''Script runs UPARSE OTU clustering. 
    Requires USEARCH by Robert C. Edgar: http://drive5.com/usearch''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
parser.add_argument('-m','--minsize', default='2', help='Min size to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--mock', default="False", help='Spike-in control: barcode label')
parser.add_argument('-mc', default='ufits_mock3.fa', help='Multi-Fasta Mock Community')
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
print "\nCMD: Quality Filtering\n%s -fastq_filter %s -fastq_maxee %s -fastqout %s\n" % (usearch, args.FASTQ, args.maxee, filter_out)
subprocess.call([usearch, '-fastq_filter', args.FASTQ, '-fastq_maxee', args.maxee, '-fastqout', filter_out], stdout = log_file, stderr = log_file)

#now run usearch8 full length dereplication
derep_out = args.out + '.EE' + args.maxee + '.derep.fa'
print "CMD: De-replication\n%s -derep_fulllength %s -sizeout -fastaout %s\n" % (usearch, filter_out, derep_out)
subprocess.call([usearch, '-derep_fulllength', filter_out, '-sizeout', '-fastaout', derep_out], stdout = log_file, stderr = log_file)

#now run usearch 8 sort by size
sort_out = args.out + '.EE' + args.maxee + '.sort.fa'
print "CMD: Sorting by Size\n%s -sortbysize %s -minsize %s -fastaout %s\n" % (usearch, derep_out, args.minsize, sort_out)
subprocess.call([usearch, '-sortbysize', derep_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = log_file, stderr = log_file)

#now run clustering algorithm
radius = str(100 - int(args.pct_otu))
otu_out = args.out + '.EE' + args.maxee + '.otus.fa'
print "CMD: Clustering OTUs\n%s -cluster_otus %s -sizein -sizeout -relabel OTU_ -otu_radius_pct %s -otus %s\n" % (usearch, sort_out, radius, otu_out)
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
    print "CMD: Chimera Filtering\n%s -uchime_ref %s -strand plus -db %s -nonchimeras %s\n" % (usearch, otu_out, uchime_db, uchime_out)
    subprocess.call([usearch, '-uchime_ref', otu_out, '-strand', 'plus', '-db', uchime_db, '-nonchimeras', uchime_out], stdout = log_file, stderr = log_file)

#option if using mock community to map OTUs to mock and add to header
if args.mock != "False":
    if args.mc == "ufits_mock3.fa":
        mock = script_path + "/lib/ufits_mock3.fa"
    else:
        mock = args.mc
    mock_file = open(mock, "r")
    mock_ref_count = 0
    for line in mock_file:
        if line.startswith (">"):
            mock_ref_count = mock_ref_count + 1
    mock_file.close()
    mock_out = args.out + '.mockmap.uc'
    print "CMD: Mapping Mock Community\n%s -usearch_global %s -strand plus -id 0.97 -db %s -uc %s\n" % (usearch, uchime_out, mock, mock_out)
    subprocess.call([usearch, '-usearch_global', uchime_out, '-strand', 'plus', '-id', '0.97', '-db', mock, '-uc', mock_out], stdout = log_file, stderr = log_file)
    #generate 2 column file from uc file
    map_out_name = args.out + '.headers.txt'
    map_out = open(map_out_name, 'ab')
    map = csv.reader(open(mock_out), delimiter='\t')
    for col in map:
        if col[-1] != "*":
            map_out.write("%s\t%s%s;pident=%s;\n" % (col[-2], col[-2], col[-1], col[3]))
    map_out.close()
    annotation = open(map_out_name, "rb")
    annotate_dict = {}
    for line in annotation:
        line = line.split("\t")
        if line:
            annotate_dict[line[0]]=line[1:]
        else:
            continue
    otu_new = args.out + '.EE' + args.maxee + '.mock.otus.fa'
    otu_update = open(otu_new, "w")
    with open(uchime_out, "r") as myfasta:
        for line in myfasta:
            if line.startswith (">"):
                line = line[1:]
                line = line.split()
                if line[0] in annotate_dict:
                    new_line = ">" + "".join(annotate_dict[line[0]])
                    otu_update.write (new_line)
                else:
                    otu_update.write (">"+ "".join(line) + "\n")
            else:
                line = re.sub('N', '', line)
                otu_update.write(line)
    otu_update.close()
    uchime_out = otu_new
    
#now map reads back to OTUs
uc_out = args.out + '.EE' + args.maxee + '.mapping.uc'
if args.map_unfiltered:
    reads = args.FASTQ
else:
    reads = filter_out
print "CMD: Mapping Reads to OTUs\n%s -usearch_global %s -strand plus -id 0.97 -db %s -uc %s\n" % (usearch, reads, uchime_out, uc_out)
subprocess.call([usearch, '-usearch_global', reads, '-strand', 'plus', '-id', '0.97', '-db', uchime_out, '-uc', uc_out], stdout = log_file, stderr = log_file)
#Build OTU table
otu_table = args.out + '.EE' + args.maxee + '.otu_table.txt'
uc2tab = script_path + "/lib/uc2otutab.py"
print "CMD: Creating OTU Table\npython %s %s > %s" % (uc2tab, uc_out, otu_table)
os.system('%s %s %s %s %s' % ('python', uc2tab, uc_out, '>', otu_table))

#run some stats on mock community if --mock option passed.
if args.mock != "False":
    KEEP_COLUMNS = ('OTUId', args.mock)
    f = csv.reader(open(otu_table), delimiter='\t')
    headers = None
    results = []
    for row in f:
        if not headers:
            headers = []
            for i, col in enumerate(row):
                if col in KEEP_COLUMNS:
                    headers.append(i)
        else:
            results.append(tuple([row[i] for i in headers]))
    num_otus = 0
    mock_found = 0
    bad_otu = []
    good_otu = []
    for row in results:
        if int(row[1]) > 0:
            num_otus = num_otus + 1
            if "pident" in row[0]:
                mock_found = mock_found + 1
                good_otu.append(int(row[1]))
            if not "pident" in row[0]:
                bad_otu.append(int(row[1]))
    spurious = num_otus - mock_found
    print "\nSummarizing data for %s" % (args.out)
    print "------------------------------------------"
    print "Total OTUs in Mock:  %i" % (mock_ref_count)
    print "Total OTUs in %s:  %i" % (args.mock, num_otus)
    print "\nReal Mock OTUs found:  %i" % (mock_found)
    good_otu = sorted(good_otu, key=int)
    print "Range of counts from Real OTUs:  %i - %i" % (good_otu[-1], good_otu[0])
    print "Lowest counts from Real OTUs:  %i, %i, %i" % (good_otu[0], good_otu[1], good_otu[2])
    print "\nSpurious OTUs found:  %i" % (spurious)
    if spurious != 0:
        bad_otu = sorted(bad_otu, key=int, reverse=True)
        print "Range of counts from Spurious OTUs:  %i - %i" % (bad_otu[0], bad_otu[-1])
        if spurious >= 3:
            print "Highest counts from Spurious OTUs:  %i, %i, %i" % (bad_otu[0], bad_otu[1], bad_otu[2])
        if spurious == 2:
            print "Highest counts from Spurious OTUs:  %i, %i" % (bad_otu[0], bad_otu[1])
        if spurious == 1:
            print "Highest count from Spurious OTUs:  %i" % (bad_otu[0])
    os.remove(map_out_name)
    os.remove(mock_out)

#Print location of files to STDOUT
print "\n------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "------------------------------------------------"
print ("Input FASTQ:           %s" % (args.FASTQ))
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

 