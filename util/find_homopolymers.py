#!/usr/bin/env python

import re
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

parser=argparse.ArgumentParser(
    description='''GC content and homopolymer length''', prog="find_homopolymers.py", usage="%(prog)s [options] -i file.fasta",
    epilog="""Written by Jon Palmer (2015) palmer.jona@gmail.com""",
    formatter_class=MyFormatter)
parser.add_argument('-i', '--fasta', dest='fasta', required=True, help='FASTA file')
parser.add_argument('-o','--out', dest='out', default='out', help='output file')
parser.add_argument('-n','--num_homopolymers', dest='homo', default='6', help='Number of homopolymers')
args=parser.parse_args()

#make default output from splitting input fasta name
if args.out == "out":
    out_base = re.split(r'\.fa', args.fasta)
    out_name = out_base[0] + ".homopolymers_" + args.homo + ".txt"
else:
    out_name = args.out

out_file = open(out_name, "w")

#print header
header_name = "SeqID"
out_file.write ("SeqID\tLength (bp)\tGC Content (%)\tHomopolymers (Len(nuc):start-stop)\n")
record_count = 0
none_count = 0
for record in SeqIO.parse(open(args.fasta, "rU"), "fasta"):
    record_count = record_count + 1
    GC_calc = GC(record.seq)
    lis = []
    for x in record.seq:
        if len(lis) != 0:
            if lis[-1][0] == x:
                lis[-1][1] += 1
            else:
                lis.append([x, 1])
        else:
            lis.append([x, 1])
    #now try to get index of those greater than argparse num_homopolymer
    num_homo = int(args.homo)
    count = 0  #keep track of index as you go through the list
    homo_out = ''
    for item in lis:
        count = count + 1
        temp_lis = []
        if item[1] >= num_homo:
            #get index from the count
            count_index = count - 1
            temp_list = lis[:count_index]
            start = sum([item[1] for item in temp_list])
            stop = start + lis[count_index][1]
            nucleotide = str(lis[count_index][0])
            if nucleotide == 'N':
                continue
            length = str(lis[count_index][1])
            homopolymers = "%s(%s):%i-%i;" % (length, nucleotide, start, stop)
            #now append to the variable
            homo_out = homo_out + homopolymers
    if homo_out == '':
        homo_out = "None found"
        none_count = none_count + 1
    else:
        homo_out = re.sub(";", "; ", homo_out)
        homo_out = re.sub("; $", ";", homo_out)
    out_file.write("%s\t%s\t%.2f\t%s\n" % (record.id, len(record.seq), GC_calc, homo_out))

homo_count = record_count - none_count
homo_pct = float(homo_count) / float(record_count) * 100
pct = "%"
print "Input Sequences: %i" % record_count
print "Seqs with homopolymers > %s: %i (%.02f%%)" % (args.homo, homo_count, homo_pct)
print "Results located here: %s" % out_name
out_file.close()