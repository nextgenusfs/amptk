#!/usr/bin/env python

#This script filters results from ufits-OTU_cluster.py
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import inspect
import subprocess
import csv
import re
from Bio import SeqIO

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-mock_filter.py', usage="%(prog)s [options] -i otu_table.txt -b BC_27\n%(prog)s -h for help menu",
    description='''Script inspects output of ufits-OTU_cluster.py and 
    determines useful threshold for OTU output based on a spike-in 
    mock community.''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--otu_table', required=True, help='Input OTU table')
parser.add_argument('-b','--mock_barcode', required=True, help='Barocde of Mock community')
parser.add_argument('-mc', dest="mock_community",default='ufits_mock3.fa', help='Multi-FASTA mock community')
parser.add_argument('--trim_data', action="store_true", help='Threshold Trim Data')

args=parser.parse_args()


def try_int(x):
    try:
        return int(x)
    except ValueError:
        return x
        
def try_subtract(x,y):
    try:
        z = x - y
        if z < 0:
            z = 0
        else:
            z = z
        return z
    except TypeError:
        return x

#check if otu_table is empty
check = os.stat(args.otu_table).st_size
if check == 0:
    print "Input file is empty"
    os.exit(1)

#get base name of files
base = re.sub('.otu_table.txt', '', args.otu_table)

#get default mock community value
if args.mock_community == "ufits_mock3.fa":
    mock = script_path + "/lib/ufits_mock3.fa"
else:
    mock = args.mock_community
#open mock community fasta and count records
mock_file = open(mock, "r")
mock_ref_count = 0
for line in mock_file:
    if line.startswith (">"):
        mock_ref_count = mock_ref_count + 1
mock_file.close()

#load in OTU table, get only OTU column and mock
KEEP_COLUMNS = ('OTUId', args.mock_barcode)
f = csv.reader(open(args.otu_table), delimiter='\t')
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
print "\nSummarizing data for %s" % (base)
print "------------------------------------------"
print "Total OTUs in Mock:  %i" % (mock_ref_count)
print "Total OTUs in %s:  %i" % (args.mock_barcode, num_otus)
print "\nReal Mock OTUs found:  %i" % (mock_found)
good_otu = sorted(good_otu, key=int)
try:
    print "Range of counts from Real OTUs:  %i - %i" % (good_otu[-1], good_otu[0])
except IndexError:
    print "\nThere does not appear to be Mock OTUs mapped in this table, run `ufits-OTU_cluster.py with the --mock parameter to generate a compatible OTU table.\n"
    os._exit(1)
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

if args.trim_data:
        threshold = raw_input("\nEnter threshold value to trim data:  ")
        num = int(threshold)
        new_table = []
        sub_table = []
        keys = []
        trim_table = []
        line_count = 0
        out_name = base + '.filtered_' + threshold + '.otu_table.txt'
        file_out = open(out_name, "wb")
        f2 = csv.reader(open(args.otu_table), delimiter='\t')
        for line in f2:
            line_count = line_count + 1
            new_table.append([try_int(x) for x in line]) #convert to integers
        for line in new_table:
            sub_table.append([try_subtract(x,num) for x in line]) #subtract threshold
        for line in sub_table:
            if max(line[1:-1]) >= 1:
                trim_table.append(line) #get rid of OTUs with only zeros
                keys.append(line[0])
        writer = csv.writer(file_out, delimiter='\t')
        writer.writerows(trim_table)
        file_out.close()
        #now lets write an updated OTU fasta file
        fasta_in = base + '.mock.otus.fa'
        fasta_out = base + '.filtered_' + threshold + '.otus.fa'
        seqs_seen = []
        count = 0
        for record in SeqIO.parse(open(fasta_in, "rU"), "fasta"):
            if record.id in keys:
                count = count + 1
                seqs_seen.append(record)
        fasta_update = open(fasta_out, "w")
        SeqIO.write(seqs_seen, fasta_update, "fasta")
        fasta_update.close()
        print "\nOTU table has been filtered to %i:\nOriginal OTUs: %i\nFiltered OTUs: %i" % (num, line_count, count)