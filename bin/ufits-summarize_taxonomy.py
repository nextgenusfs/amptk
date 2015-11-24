#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import csv, argparse, re, os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from lib.stackedBarGraph import StackedBarGrapher as StackedBarGrapher

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='ufits-summarize_taxonomy.py', usage="%(prog)s -i ufits.otu_table.taxonomy.txt\n%(prog)s -h for help menu",
    description='''Script that produces summary figures and data tables from UFITS taxonomy OTU table.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--table', dest="table", required=True, help='OTU Table (Required)')
parser.add_argument('-o','--out', dest="out", default='ufits-summary', help='Base name of Output Files')
parser.add_argument('--graphs', dest="graphs", action="store_true", help='Create Stacked Bar Graphs')
parser.add_argument('--format', dest="format", default='eps', choices=['eps','svg','png','pdf'], help='Image format')
parser.add_argument('--percent', dest="percent", action="store_true", help='Convert to Pct of Sample')
parser.add_argument('--font_size', dest="size", default=8, help='Font Size')
args=parser.parse_args()

def try_int(x):
    try:
        return int(x)
    except ValueError:
        return x       

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

def tally(x):
    from collections import Counter
    tally = Counter()
    for i in x:
        tally[i] += 1
    return tally

def getClassCounts(table):
    global k,p,c,o,f,g,s
    k = []
    p = []
    c = []
    o = []
    f = []
    g = []
    s = []   
    for item in table:
        tax = item
        tax = tax.replace(';', ',')
        tax = tax.split(",")
        for i in tax:
            if i.startswith("k:"):
                k.append(i.rsplit(":", 1)[-1].split(" (")[0])
            if i.startswith("p:"):
                p.append(i.rsplit(":", 1)[-1].split(" (")[0])
            if i.startswith("c:"):
                c.append(i.rsplit(":", 1)[-1].split(" (")[0])
            if i.startswith("o:"):
                o.append(i.rsplit(":", 1)[-1].split(" (")[0])
            if i.startswith("f:"):
                f.append(i.rsplit(":", 1)[-1].split(" (")[0])
            if i.startswith("g:"):
                g.append(i.rsplit(":", 1)[-1].split(" (")[0])
            if i.startswith("s:"):
                s.append(i.rsplit(":", 1)[-1].split(" (")[0])
    return k,p,c,o,f,g,s

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

def processTax(uniq, L, name):
    FinalTab = []
    OutputTable = []
    for i in range(len(allSamples)):
        finalLine = []
        printLine = [allSamples[i]]
        for x in range(len(uniq)):
            countTax = L[i].get(uniq[x]) or 0
            finalLine.append(countTax)
            printLine.append(countTax)
        FinalTab.append(finalLine)
        OutputTable.append(printLine)
    Tax = ['Sample']
    Tax.append(uniq)
    Tax = flatten(Tax)
    OutputTable.insert(0,Tax)
    table_out = args.out + '.tax.' + name + '.csv'
    with open(table_out, 'w') as tableOut:
        for line in OutputTable:
            new_line = ','.join(str(x) for x in line)+'\n'
            tableOut.write(new_line)    
    #check if not drawing charts
    if args.graphs:
        percent_table = []
        if args.percent:
            for line in FinalTab:
                new_line = []
                sums = sum(line)
                new_line.append([x/float(sums) * 100 for x in line])
                new_line = flatten(new_line)
                percent_table.append(new_line)
        else:
            percent_table = FinalTab
        out = args.out + '.tax.' + name + '.' + args.format
        SBG = StackedBarGrapher()
        #load pd df
        df = pd.DataFrame(percent_table)
        #work on colors, 12 preferred, add to those if necessary
        pref_colors=["#023FA5","#7D87B9","#BEC1D4","#D6BCC0",
                "#BB7784","#4A6FE3","#8595E1","#B5BBE3",
                "#E6AFB9","#E07B91","#D33F6A","#11C638","#8DD593",
                "#C6DEC7","#EAD3C6","#F0B98D","#EF9708","#0FCFC0",
                "#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4"]
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
        d_labels = allSamples
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if args.percent:
            YLabel = "Percent of Community"
        else:
            YLabel = "Number of OTUs per taxonomy level"
        SBG.stackedBarPlot(ax,df,d_colors,edgeCols=['#000000']*length,xLabels=allSamples,gap=0.25,endGaps=True,xlabel="Samples", ylabel=YLabel)
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
        ax.get_xticklabels().set_fontsize(int(args.font_size))

        #setup the plot
        fig.subplots_adjust(bottom=0.4)
        fig.set_tight_layout(True) 
        fig.savefig(out, format=args.format, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)



sub_table = []
f2 = csv.reader(open(args.table), delimiter='\t')
for line in f2:
    sub_table.append([try_int(x) for x in line]) #convert to integers
    
header = sub_table[:1] #pull out header row
header = header[0] #just get first list, double check
headerCount = len(header) - 1

allSamples, Lk, Lp, Lc, Lo, Lf, Lg, Ls, Tk, Tp, Tc, To, Tf, Tg, Ts = ([] for i in range(15))
for i in range(1,headerCount):
    allSamples.append(header[i])
    tab = []
    OTUs = 0
    for line in sub_table:
        if line[i] > 0:
            OTUs += 1
            tab.append(line[-1])
    
    getClassCounts(tab)    
    Kingdom = tally(k)
    Kingdom['unclassified'] = OTUs - sum(Kingdom.values())
    Lk.append(dict(Kingdom))
    Tk.append(list(set(Kingdom)))    
    Phylum = tally(p)
    Phylum['unclassified'] = OTUs - sum(Phylum.values())
    Lp.append(dict(Phylum))
    Tp.append(list(set(Phylum)))
    Class = tally(c)
    Class['unclassified'] = OTUs - sum(Class.values())
    Lc.append(dict(Class))
    Tc.append(list(set(Class)))
    Order = tally(o)
    Order['unclassified'] = OTUs - sum(Order.values())
    Lo.append(dict(Order))
    To.append(list(set(Order)))
    Family = tally(f)
    Family['unclassified'] = OTUs - sum(Family.values())
    Lf.append(dict(Family))
    Tf.append(list(set(Family)))
    Genus = tally(g)
    Genus['unclassified'] = OTUs - sum(Genus.values())
    Lg.append(dict(Genus))
    Tg.append(list(set(Genus)))   
    '''
    Species = tally(s)
    Species['unclassified'] = OTUs - sum(Species.values())
    Ls.append(dict(Species))
    Ts.append(list(set(Species)))
    '''
#kingdom
uniqK = list(set(flatten(Tk)))
uniqK.sort()
processTax(uniqK, Lk, 'kingdom')

#phylum
uniqP = list(set(flatten(Tp)))
uniqP.sort()
processTax(uniqP, Lp, 'phylum')

#class
uniqC = list(set(flatten(Tc)))
uniqC.sort()
processTax(uniqC, Lc, 'class')

#order
uniqO = list(set(flatten(To)))
uniqO.sort()
processTax(uniqO, Lo, 'order')

#family
uniqF = list(set(flatten(Tf)))
uniqF.sort()
processTax(uniqF, Lf, 'family')

#genus
uniqG = list(set(flatten(Tg)))
uniqG.sort()
processTax(uniqG, Lg, 'genus')

'''#species
uniqS = list(set(flatten(Ts)))
uniqS.sort()
processTax(uniqS, Ls, 'species')
'''