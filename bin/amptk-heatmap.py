#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import argparse
import os
import sys
from natsort import natsorted
try:
    import seaborn as sns
except:
    pass

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='amptk-heatmap.py', usage="%(prog)s -i amptk.otu_table.txt\n%(prog)s -h for help menu",
    description='''Draw heatmap from OTU table generated with AMPtk.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

map_colors = ['Accent','BrBG','BuGn','Dark2','PRGn','Pastel1', 'Pastel2','Spectral','bwr','copper','cubehelix','gist_yarg','gist_earth_r','gist_heat','Blues','Reds','Purples','Greens','BuPu','BuGn','OrRd','YlGnBu','YlOrRd','terrain', 'brg', 'gnuplot','bone_r','terrain_r','sprint','summer','ocean_r']

parser.add_argument('-i','--table', dest="table", required=True, help='OTU Table (Required)')
parser.add_argument('-o','--output', dest="output", required=True, help='Output File (Required)')
parser.add_argument('--format', dest="format", default='eps', choices=['eps','svg','png','pdf'], help='Image format')
parser.add_argument('--col_order', dest="col_order", default="naturally", help='Provide comma separated list [naturally, None, list]')
parser.add_argument('--font_size', dest="size", default=10, help='Font Size')
parser.add_argument('--border_color', dest="border_color", default='lightgray', choices=['black', 'white','lightgray'], help='color of border between cells')
parser.add_argument('--zero_color', dest="zero_color", default='white', choices=['white','lightgray','black','snow'], help='Color for zeros')
parser.add_argument('--square', dest="square", action="store_true", help='Maintain aspect ratio')
parser.add_argument('--colors', dest="colors", default='YlOrRd', choices=map_colors, help='Choose color palette')
parser.add_argument('--percent', dest="percent", action="store_true", help='Convert to Pct of Sample')
parser.add_argument('--min_num', dest="min_num", default=1, help='Min number of reads per OTU')
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

def try_int(x):
    try:
        return int(x)
    except ValueError:
        return x

def convert_binary(x):
    if isinstance(x, float):
        if x >= 0.3:
            return 1
        else:
            return 0
    else:
        return x

if args.table.rsplit(".", 1)[1] == 'txt':
    delim = "\t"
if args.table.rsplit(".", 1)[1] == 'csv':
    delim = ","
    
#open OTU table
#get the OTU header info (depending on how OTU table was constructed, this might be different, so find it as you need for indexing)
with open(args.table, 'rU') as f:
    first_line = f.readline()
    OTUhead = first_line.split(delim)[0]
    if 'Taxonomy' in first_line.split(delim)[-1]:
        tax = True
    else:
        tax = False

#this is so dumb, just use pandas
df = pd.read_csv(args.table, sep=delim, index_col=0)
if 'Taxonomy' in df.columns.values:
    del df['Taxonomy']
print(df)
sys.exit(1)

f2 = csv.reader(open(args.table), delimiter=delim)
taxTable = []
#test for taxonomy
if tax:
    for line in f2:
        del line[-1]
        taxTable.append(line)
else:
    for line in f2:
        taxTable.append(line)
print(taxTable)
otuDict = {rows[0]:rows[-1] for rows in taxTable}

#first convert the numbers to integers instead of strings
new_table = []
temp_table = []
min_table = []
binary_table = []
for line in taxTable:
    new_table.append([try_int(x) for x in line]) #convert to integers
header = new_table[:1]
header = header[0]
    
if args.min_num > 1:
    min_table.append(header)
    for line in new_table[1:]:
        if sum(line[1:-1]) >= int(args.min_num):
            min_table.append(line)
else:
    min_table = new_table

if args.percent:
    #get sums of columns, insert into second row of table
    for line in min_table[1:]:
        temp_table.append(line[1:])
    sums = [sum(i) for i in zip(*temp_table)]
    sums.insert(0,OTUhead)

    #now get percentage for each
    for line in min_table[1:]:
        new_line = []
        new_line.insert(0,line[0])
        new_line.append([x/y * 100 for x,y in zip(line[1:], sums[1:])])     
        new_line = flatten(new_line)
        binary_table.append(new_line)
    binary_table.insert(0,header)
else:
    binary_table = min_table

#Now sort the table by row name naturally
if args.col_order == "naturally":
    sortHead = natsorted(header)  #now sort the header row
elif args.col_order == "None":
    sortHead = header
else:
    sortHead = args.col_order.split(",")
    sortHead.append(OTUhead)
otu = sortHead.index(OTUhead) #pull out the OTU header
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

if count >= 100:
    print("Your table has %i rows (OTUs)" % count)
    print("It is likely not going to be displayed very well")
    filter_num = input('How many rows should I print: ')
    filter_num = int(filter_num) + 1
    finalTable = finalTable[:filter_num]

#Ok, now convert the sorted table to a panda data frame
Cols = finalTable[0][1:] #get column names
Index = [] 
Index.append([item[0] for item in finalTable]) #remove first column
del Index[0][0]
finalTable = finalTable[1:]
for x in finalTable:
    del x[0]

for item in Index[0]:
    taxon = otuDict.get(item)
    #print item,taxon.split(",")[-1] #print the OTU and lowest taxonomic classification
    
#construct panda dataframe with appropriate headers and index
df = pd.DataFrame(finalTable, index=Index, columns=Cols)
print(df)
#get sizes
width = len(df.columns)/4
height = len(df.index)/4

# Plot it out
#get color palette from argparse
cmap = plt.get_cmap(args.colors)
cmap.set_under(args.zero_color)
fig, ax = plt.subplots(figsize=(width,height))
heatmap = ax.pcolor(df, cmap=cmap, vmin=0.0001, linewidths=1, edgecolors=args.border_color)

cbar = plt.colorbar(heatmap, cmap=cmap, cax=None, ax=None, shrink=0.4)
if args.percent:
    cbar.ax.set_ylabel('Percent of reads per sample')
else:
    cbar.ax.set_ylabel('# of reads')


# Format
fig = plt.gcf()
fig.set_size_inches(8,11)

# turn off the frame
ax.set_frame_on(False)

# put the major ticks at the middle of each cell
ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

# Set the labels
ax.set_xticklabels(df.columns, minor=False)
ax.set_yticklabels(df.index, minor=False)

# rotate the
plt.xticks(rotation=90)
ax.grid(False)

# Turn off all the ticks
ax = plt.gca()

for t in ax.xaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False
for t in ax.yaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False

#set the font size - i wish I knew how to do this proportionately.....
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(args.size)

if args.square:
    #make axis equal
    ax.set_aspect('equal')

#plt.show()
plt.savefig(args.output, format=args.format, dpi=1000)
