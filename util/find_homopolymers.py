#!/usr/bin/env python

import re
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

parser=argparse.ArgumentParser(
    description='''GC content and homopolymer length''', prog="find_homopolymers.py", usage="%(prog)s [options] [options] file.fasta",
    epilog="""Written by Jon Palmer (2015) palmer.jona@gmail.com""",
    formatter_class=MyFormatter)
parser.add_argument('fasta', help='FASTA file')
parser.add_argument('-n','--num_homopolymers', dest='homo', default='6', help='Number of homopolymers')
parser.add_argument('--print_pretty', action='store_true', help='print to terminal in aligned columns')
args=parser.parse_args()

#print header
header_name = "SeqID"
if args.print_pretty:
    print header_name.ljust(20) + "\tLength (bp)\tGC Content (%)\tHomopolymers (Len(nuc):start-stop)"
else:
    print "SeqID\tLength (bp)\tGC Content (%)\tHomopolymers (Len(nuc):start-stop)"

for record in SeqIO.parse(open(args.fasta, "rU"), "fasta"):
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
            length = str(lis[count_index][1])
            homopolymers = "%s(%s):%i-%i;" % (length, nucleotide, start, stop)
            #now append to the variable
            homo_out = homo_out + homopolymers
    if homo_out == '':
        homo_out = "None found"
    else:
        homo_out = re.sub(";", "; ", homo_out)
        homo_out = re.sub("; $", ";", homo_out)
    if args.print_pretty:
        print "%-20s\t%-8s\t%-10.2f\t%s" % (record.id, len(record.seq), GC_calc, homo_out)
    else:
        print "%s\t%s\t%.2f\t%s" % (record.id, len(record.seq), GC_calc, homo_out)