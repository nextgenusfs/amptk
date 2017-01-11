#!/usr/bin/env python

import csv
import argparse
import re
import os
import sys
from natsort import natsorted

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='amptk-heatmap.py', usage="%(prog)s -i amptk.otu_table.txt\n%(prog)s -h for help menu",
    description='''''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--table', dest="table", required=True, help='OTU Table (Required)')
parser.add_argument('-p','--percent_threshold', dest="percent", default='0.3', help='OTU Threshold (percent)')
parser.add_argument('--only_convert_binary', dest="binary", action='store_true', help='Only convert to binary')
parser.add_argument('--col_order', dest="col_order", default="naturally", help='Provide comma separated list')
args=parser.parse_args()

percent = float(args.percent)

def try_int(x):
    try:
        return int(x)
    except ValueError:
        return x

def convert_binary(x):
    if isinstance(x, float):
        if x >= percent:
            return 1
        else:
            return 0
    else:
        return x
def only_convert_binary(x):
    if isinstance(x, int):
        if x >= 1:
            return 1
        else:
            return 0
    else:
        return x

def convert_percent(x,y):
    try:
        percent = float(x) / float(y)
        percent = percent * 100
        return percent
    except ValueError:
        return x
        
#open OTU table   
f2 = csv.reader(open(args.table), delimiter='\t')

if args.binary:
    #first convert the numbers to integers, then to presence/absence
    new_table = []
    temp_table = []
    binary_table = []
    for line in f2:
        new_table.append([try_int(x) for x in line]) #convert to integers
    #print new_table
    for line in new_table:
        if line[0] == 'OTUId':
            binary_table.append(line)
        else:
            binary_table.append([only_convert_binary(x) for x in line]) #convert to binary

else:
    #first convert the numbers to integers, then to presence/absence
    new_table = []
    temp_table = []
    binary_table = []
    for line in f2:
        new_table.append([try_int(x) for x in line]) #convert to integers
    #print new_table
    for line in new_table:
        if line[0] == 'OTUId':
            temp_table.append(line)
        else:
            subline = line[1:]
            OTU_sum = sum(subline)
            temp_table.append([convert_percent(x,OTU_sum) for x in line]) #convert to percent
    for line in temp_table:
        binary_table.append([convert_binary(x) for x in line]) #convert to binary

#Now sort the table by row name naturally
header = binary_table[:1] #pull out header row
header = header[0] #just get first list, double check
if args.col_order == "naturally":
    sortHead= natsorted(header)  #now sort the header row
else:
    sortHead = args.col_order.split(",")
    sortHead.append('OTUId')
otu = sortHead.index('OTUId') #pull out the OTU header
sortHead.insert(0, sortHead.pop(otu)) #re-insert the OTU header in first column

#sort the table without header row and then add back
sortTable = natsorted(binary_table[1:]) #sort the table without header row
sortTable.insert(0,header) #now add back header row

#get index of the sorted header list 
listIndex = []
for i in sortHead:
    listIndex.append(header.index(i))

#finally loop through the sorted table, and re-order the data
finalTable = []
count = 0
for line in sortTable:
    count += 1
    lineList= []
    for item in listIndex:
        lineList.append(line[item])
    finalTable.append(lineList)

for line in finalTable:
    print ",".join(str(x) for x in line)