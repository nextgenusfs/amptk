#!/usr/bin/env python

#script to pivot OTU table and append to metadata

import sys, os, itertools, csv, argparse
from natsort import natsorted

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-remove_samples.py',
    description='''Script to sub-sample reads down to the same number for each sample (barcode)''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-m','--meta', required=True, help='Meta data (csv format, e.g. from excel)')
parser.add_argument('-i','--input', required=True, help='OTU table')
parser.add_argument('-o','--out', required=True, help='Output name')
args=parser.parse_args()


def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList 
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

#first import OTU table, pivot and convert to dictionary
with open(args.input,'rU') as infile:
    if args.input.endswith('.csv'):
        d = ','
    else:
        d = '\t'
    reader = itertools.izip(*csv.reader(infile, delimiter=d))
    otuDict = {rows[0]:rows[1:] for rows in reader}

#print otuDict
#now load meta table and append tuple from dictionary
finalTable = []
with open(args.meta, 'rU') as input:
    reader = csv.reader(input)
    count = 0
    for line in reader:
        count += 1
        if count == 1:
            header = otuDict.get('OTUId')
            header = list(header)
            line.append(header)
        else:
            otus = otuDict.get(line[0].rstrip()) or "Sample not found"
            if otus != "Sample not found":
                otus = list(otus)
            else:
                print "%s not found in OTU table, they must match exactly!" % line[0]
            line.append(otus)
        line = flatten(line)
        finalTable.append(line)

#save to file
with open(args.out, 'w') as output:
    for line in finalTable:
        if line[-1] != "Sample not found":
            join_line = (','.join(str(x) for x in line))
            output.write("%s\n" % join_line)