#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import matplotlib
matplotlib.use('agg')
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
import pandas as pd
import numpy as np
import argparse
import os
import sys
from collections import OrderedDict
from natsort import natsorted
from amptk import amptklib
from amptk import stackedBarGraph

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

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
    
def get_colors(num_colors):
    import colorsys
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        clr = colorsys.hls_to_rgb(hue, lightness, saturation)
        hex = '#%02x%02x%02x' % (int(clr[0]*256), int(clr[1]*256), int(clr[2]*256))
        hex.upper()
        colors.append(hex)
    return colors
  
def drawBarGraph(df, output, args=False):
        SBG = stackedBarGraph.StackedBarGrapher()
        #work on colors, 12 preferred, add to those if necessary
        pref_colors=["#023FA5","#7D87B9","#BEC1D4","#D6BCC0",
                "#BB7784","#4A6FE3","#8595E1","#B5BBE3",
                "#E6AFB9","#E07B91","#D33F6A","#11C638","#8DD593",
                "#C6DEC7","#EAD3C6","#F0B98D","#EF9708","#0FCFC0",
                "#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4"]
        uniq = df.columns.values.tolist()
        length = len(uniq)
        if length > len(pref_colors):
            extra = length - len(pref_colors) - 1
            if extra > 0:
                add_colors = get_colors(extra)
                pref_colors.append(add_colors)
                pref_colors.append('#ededed') 
            else:
                pref_colors.append('#ededed')
        else:
            cut = len(uniq) - 1
            pref_colors = pref_colors[:cut]
            pref_colors.append('#ededed')
        d_colors = flatten(pref_colors)

        #labels
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if args.percent:
            YLabel = "Percent of Community"
        else:
            YLabel = "Number of OTUs per taxonomy level"
        SBG.stackedBarPlot(ax,df,d_colors,edgeCols=['#000000']*length,xLabels=df.index.values.tolist(),gap=0.25,endGaps=True,xlabel="Samples", ylabel=YLabel)
        plt.title("Taxonomy Summary")
        #get the legend
        legends = [] 
        i = 0 
        for column in df.columns: 
            legends.append(mpatches.Patch(color=d_colors[i], label=uniq[i])) 
            i+=1 
        lgd = ax.legend(handles=legends, fontsize=6, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(True)
        #set the font size - i wish I knew how to do this proportionately.....
        for item in ax.get_xticklabels():
            item.set_fontsize(args.size)

        #setup the plot
        fig.subplots_adjust(bottom=0.4)
        fig.savefig(output, format=args.format, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)


def main(args):
    parser=argparse.ArgumentParser(prog='amptk-summarize_taxonomy.py', usage="%(prog)s -i amptk.otu_table.taxonomy.txt\n%(prog)s -h for help menu",
        description='''Script that produces summary figures and data tables from AMPtk taxonomy OTU table.''',
        epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('-i','--table', dest="table", required=True, help='OTU Table (Required)')
    parser.add_argument('-o','--out', dest="out", default='amptk-summary', help='Base name of Output Files')
    parser.add_argument('--graphs', dest="graphs", action="store_true", help='Create Stacked Bar Graphs')
    parser.add_argument('--format', dest="format", default='eps', choices=['eps','svg','png','pdf'], help='Image format')
    parser.add_argument('--percent', dest="percent", action="store_true", help='Convert to Pct of Sample')
    parser.add_argument('--counts', default='binary', choices=['actual', 'binary'], help='Use actual abundance or binary')
    parser.add_argument('--font_size', dest="size", type=int, default=8, help='Font Size')
    args=parser.parse_args(args)
    
    #rewrite as I don't know what is going on with old script and why it is now failing
    #goal is just to loop through the taxonomy and create CSV and optional graph of each taxonomy level
    taxLookup = {'d': 'domain', 'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o':'order', 'f': 'family', 'g': 'genus', 's': 'species'}
    
    #lets create an OTU : taxonomy dictionary and then a sample : OTU dictionary
    otuDict = OrderedDict()
    sampleDict = OrderedDict()
    delim = '\t'
    header = []
    with open(args.table, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('#OTU'): #header row
                if line.count('\t') > line.count(','):
                    delim = '\t'
                else:
                    delim = ','
                header = line.split(delim)
                if not 'Taxonomy' in header:
                    print('Error: Taxonomy not found in header of OTU table')
                    sys.exit(1)
                for i,x in enumerate(header):
                    if i == 0 or x == 'Taxonomy' or x == 'taxonomy':
                        continue
                    sampleDict[i] = {'sample': x, 'OTUs': {}}   
                continue
            cols = line.split(delim)
            #populate OTU dict
            #GS|100.0|GQ253122_S001576427;k:Bacteria,p:"Proteobacteria",c:Alphaproteobacteria,o:Sphingomonadales,f:Sphingomonadaceae,g:Sphingomonas,s:Sphingomonas_glacialis;
            if not cols[-1].startswith('No'):
                score, tax = cols[-1].split(';', 1)
                tax = tax.rstrip(';')
                otuDict[cols[0]] = {'score': score, 'domain': None, 'kingdom': None, 'phylum': None, 'class': None, 'order': None,
                                    'family': None, 'genus': None, 'species': None}
                for w in tax.split(','):
                    try:
                        level, value = w.split(':')
                        otuDict[cols[0]][taxLookup[level]] = value
                    except ValueError:
                        pass
            #now populate sample dictionary
            for z,y in enumerate(cols):
                if z == 0 or z == len(cols)-1:
                    continue
                sampleDict[z]['OTUs'][cols[0]] = int(y)
    
    #we can now loop through each taxonomic level constructing the taxonomy summary into pandas dataframe
    for taxLevel in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']: #, , 'species']:
        taxSummary = OrderedDict()
        for k,v in sampleDict.items():
            sampleID = v['sample']
            for t in natsorted(v['OTUs'].items()):
                otuName, Count = t
                if args.counts == 'binary':
                    if Count > 0:
                        Count = 1
                taxInfo = otuDict[otuName][taxLevel]
                #print(sampleID, otuName, Count, taxInfo)
                if not taxInfo:
                    taxInfo = 'Unclassified'
                if not sampleID in taxSummary:
                    taxSummary[sampleID] =  {taxInfo: Count}
                else:
                    if not taxInfo in taxSummary[sampleID]:
                        taxSummary[sampleID][taxInfo] = Count
                    else:
                        taxSummary[sampleID][taxInfo] += Count
        df = pd.DataFrame(taxSummary)
        if len(df.index.values.tolist()) == 1 and 'Unclassified' == df.index.values.tolist()[0]:
            continue
        if args.percent:
            df = (100. * df / df.sum())
        df.transpose().to_csv(args.out+'.'+taxLevel+'.csv')
        if args.graphs:
            drawBarGraph(df.transpose(), args.out+'.'+taxLevel+'.'+args.format, args=args)
            
if __name__ == "__main__":
    main(sys.argv[1:])
