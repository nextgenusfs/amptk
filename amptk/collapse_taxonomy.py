#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import pandas as pd
from natsort import natsorted

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='amptk-collapse_taxonomy.py',
    description='''.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input OTU table')
parser.add_argument('--tax_levels', default='g', choices=['k','p','c','o','f','g','s'], help='Levels of taxonomy needed to collapse')
parser.add_argument('-o','--out', required=True, help='Basename for output files')
args=parser.parse_args()

taxlevel = {'k': (0, 'kingdom'), 'p': (1, 'phylum'),'c': (2, 'class'),'o': (3, 'order'),'f': (4, 'family'),'g': (5,'genus'),'s': (6, 'species')}

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

if args.input.endswith('.csv'):
    delim = ','
else:
    delim = '\t'

df = pd.read_csv(args.input, sep=delim)
otuDict = {}
notax = {}
taxlevelfilter = taxlevel.get(args.tax_levels)[0]
taxlevelname = taxlevel.get(args.tax_levels)[1]
for i in range(0, len(df)):
    otu = df['OTUId'][i]
    tax = df['Taxonomy'][i].split(';')[-1]
    if tax.count(',') < taxlevelfilter: #need to have at least x levels of tax
        if not tax in notax:
            notax[tax] = [otu]
        continue
    if not tax in otuDict:
        otuDict[tax] = [otu]
    else:
        otuDict[tax].append(otu)

singletons = []
for k,v in list(notax.items()):
    singletons.append(v[0])

#re-index to OTUId to make retrieving easier
df.set_index('OTUId', inplace=True)
print("Parsing %i OTUs" % len(df))
new = pd.DataFrame(df, index=singletons)
uniqcount = len(new)
print("Found %i OTUs with unique taxonomy to %s" % (len(new), taxlevelname))
for k,v in list(otuDict.items()):
    if len(v) == 1:
        df2 = pd.DataFrame(df, index=v)
        new = pd.concat([new,df2])
    else:
        df3 = pd.DataFrame(df, index=v)
        df4 = df3.sum(axis=0, numeric_only=True)
        df4['Taxonomy'] = k
        otus = [x.replace('OTU_', '') for x in v]
        otu_name = 'cOTUs_'+'_'.join(otus)
        new.loc[otu_name] = df4

#sort the table
new2 = new.reindex(index=natsorted(new.index))
new3 = new2.reindex(columns=natsorted(new2.columns))

print("Collapsing %i OTUs into %i with identical taxonomy to %s levels" % (len(df)-uniqcount, len(new3) - uniqcount, taxlevelname))

#create output files
otu_out = args.out+'.csv'
new3.to_csv(otu_out, sep=',')
with open(args.out+'.collapsed_tax.txt', 'w') as output:
    output.write('%s\t%s\n' % ('Taxonomy', 'OTU(s)'))
    for k,v in list(notax.items()):
        if 'tax=' in k:
            k = k.replace('tax=', '')
        output.write('%s\t%s\n' % (k, ', '.join(v)))
    for k,v in list(otuDict.items()):
        if 'tax=' in k:
            k = k.replace('tax=', '')
        output.write('%s\t%s\n' % (k, ', '.join(v)))
        

