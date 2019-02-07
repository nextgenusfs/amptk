#!/usr/bin/env python

from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
import sys, os, itertools, csv, argparse
from natsort import natsorted

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

def countOTUs(file):
    count = 0
    with open(file, 'r') as input:
        for line in input:
            if line.startswith('#OTU ID'):
                continue
            else:
                count += 1
    return count

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

def all_indices(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices
    
def transposeTable(metadata, otu_Dict, outfile, filter):
    #now load meta table and append tuple from dictionary
    global errors
    finalTable = []
    errors = []
    with open(metadata, 'r') as input:
        reader = csv.reader(input)
        count = 0
        for line in reader:
            count += 1
            if count == 1:
                header = otu_Dict.get('#OTU ID')
                header = list(header)
                if filter:
                    result = []
                    for i in filter:
                        result.append(header[i])
                else:
                    result = header
                line.append(result)
            else:
                otus = otu_Dict.get(line[0].rstrip()) or "Sample not found"
                if otus != "Sample not found":
                    otus = list(otus)
                    if filter:
                        result2 = []
                        for i in filter:
                            result2.append(otus[i])
                    else:
                        result2 = otus
                    line.append(result2)
                else:
                    line.append(otus)
                    errors.append(line[0])
            line = flatten(line)
            finalTable.append(line)
    #save to file
    with open(outfile, 'w') as output:
        for line in finalTable:
            if line[-1] != "Sample not found":
                join_line = (','.join(str(x) for x in line))
                output.write("%s\n" % join_line)

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-merge_metadata.py',
		description='''Takes a meta data csv file and OTU table and makes transposed output files.''',
		epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)
	parser.add_argument('-m','--meta', required=True, help='Meta data (csv format, e.g. from excel)')
	parser.add_argument('--split_taxonomy', choices=['k','p','c','o','f','g'], help='Split output files based on taxonomy levels')
	parser.add_argument('-i','--input', required=True, help='OTU table')
	parser.add_argument('-o','--out', required=True, help='Output name')
	args=parser.parse_args(args)

	if not args.out.endswith('.csv'):
		print("Error: output file must end with .csv")
		sys.exit(1)

	#first import OTU table, pivot and convert to dictionary
	with open(args.input,'r') as infile:
		if args.input.endswith('.csv'):
			d = ','
		else:
			d = '\t'
		reader = zip(*csv.reader(infile, delimiter=d))
		otuDict = {rows[0]:rows[1:] for rows in reader}

	taxonomy = []
	if args.split_taxonomy:
		#get taxonomy options into list
		try:
			for x in otuDict['Taxonomy']:
				try:
					taxlv = x.split(args.split_taxonomy+':')[1]
					taxlv = taxlv.split(',')[0]
					taxlv = taxlv.split(' ')[0]
				except IndexError:
					continue
				if not taxlv in taxonomy:
					taxonomy.append(taxlv)
		except KeyError:
			print("Error: OTU table does not contain taxonomy information")
			sys.exit(1)
		#now have list of search terms to use in the taxonomy list to filter the table
		#print taxonomy

	#run entire dataset first, then if taxonomy given, loop through each filtering via index lookup
	num_otus = countOTUs(args.input)
	print("Working on complete OTU table: found %i OTUs" % num_otus)
	transposeTable(args.meta, otuDict, args.out, False)
	if taxonomy:
		for y in taxonomy:
			index_match = [i for i, j in enumerate(otuDict['Taxonomy']) if args.split_taxonomy+':'+y in j]
			if '/' in y:
				y = y.replace('/', '_')
			tmpout = args.out.split('.csv')[0]+'.'+y+'.csv'
			print("Working on table for %s: found %i OTUs" % (y, len(index_match)))
			transposeTable(args.meta, otuDict, tmpout, index_match)
		
	if len(errors) > 0:
		print("ERROR %s: not found in OTU table. Names must max exactly." % ', '.join(errors))
		sys.exit(1)
		
	print(len(otuDict['Taxonomy']))
	for i in range(0,len(otuDict['Taxonomy'])):
		if 'Basidiomycota' in otuDict['Taxonomy'][i]:
			print(i, otuDict['OTUId'][i], otuDict['Taxonomy'][i])

if __name__ == "__main__":
	main(args)