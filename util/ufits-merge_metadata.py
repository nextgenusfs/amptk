#!/usr/bin/env python

import sys, os, itertools, csv
from natsort import natsorted

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


#script to pivot OTU table and append to metadata

if len(sys.argv) < 2:
    print "Usage: %s metadata.csv OTU_table output" % (sys.argv[0])
    sys.exit(1)

#first import OTU table, pivot and convert to dictionary of tuples?

with open(sys.argv[2],'rU') as infile:
    reader = itertools.izip(*csv.reader(infile, delimiter='\t'))
    otuDict = {rows[0]:rows[1:] for rows in reader}

#print otuDict

#now load meta table and append tuple from dictionary
finalTable = []
with open(sys.argv[1], 'rU') as input:
    reader = csv.reader(input)
    count = 0
    for line in reader:
        count += 1
        if count == 1:
            header = otuDict.get('OTUId')
            header = list(header)
            line.append(header)
        else:
            otus = otuDict.get(line[0]) or "Sample not found"
            if otus != "Sample not found":
                otus = list(otus)
            line.append(otus)
        line = flatten(line)
        finalTable.append(line)

#save to file
with open(sys.argv[3], 'w') as output:
    for line in finalTable:
        if line[-1] != "Sample not found":
            join_line = (','.join(str(x) for x in line))
            output.write("%s\n" % join_line)