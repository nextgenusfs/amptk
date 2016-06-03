#!/usr/bin/env python

#This script filters results from ufits-OTU_cluster.py
#written by Jon Palmer palmer.jona at gmail dot com

import os, argparse, inspect, subprocess, csv, logging, sys, math
from Bio import SeqIO
from natsort import natsorted
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.ufitslib as ufitslib


class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-filter.py',
    description='''Script inspects output of ufits-OTU_cluster.py and 
    determines useful threshold for OTU output based on a spike-in 
    mock community.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--otu_table', required=True, help='Input OTU table')
parser.add_argument('-f','--fasta', required=True, help='Input OTUs (multi-fasta)')
parser.add_argument('-b','--mock_barcode', help='Barocde of Mock community')
parser.add_argument('-p','--index_bleed',  type=float, help='Index Bleed filter. Defaul: auto')
parser.add_argument('-s','--subtract', default=0, help='Threshold to subtract')
parser.add_argument('-n','--normalize', default='y', choices=['y','n'], help='Normalize OTU table prior to filtering')
parser.add_argument('--mc',default='synmock', help='Multi-FASTA mock community')
parser.add_argument('-d','--delimiter', default='csv', choices=['csv','tsv'], help='Delimiter')
parser.add_argument('--col_order', dest="col_order", default="naturally", help='Provide comma separated list')
parser.add_argument('--keep_mock', action='store_true', help='Keep mock sample in OTU table (Default: False)')
parser.add_argument('-o','--out', help='Base output name')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

if not args.out:
    #get base name of files
    base = args.otu_table.split(".otu_table")
    base = base[0]
else:
    base = args.out

#remove logfile if exists
log_name = base + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
ufitslib.log.info("Operating system: %s, %s" % (sys.platform, ufitslib.get_version()))

#get usearch location/name
usearch = args.usearch
try:
    usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
except OSError:
    ufitslib.log.warning("%s not found in your PATH, exiting." % usearch)
    os._exit(1)
ufitslib.log.info("USEARCH version: %s" % usearch_test)

#check if otu_table is empty
ufitslib.log.info("Loading OTU table: %s" % args.otu_table)
check = os.stat(args.otu_table).st_size
if check == 0:
    ufitslib.log.error("Input OTU table is empty")
    os._exit(1)

if args.delimiter == 'csv':
    delim = ','
    ending = '.csv'
elif args.delimiter == 'tsv':
    delim = '\t'
    ending = '.txt'

#setup outputs
sorted_table = base+'.sorted'+ending
normal_table_pct = base+'.normalized.pct'+ending
normal_table_nums = base+'.normalized.num'+ending
subtract_table = base+'.normalized.subtract'+ending
final_table = base+'.final'+ending
final_binary_table = base+'.final.binary'+ending
stats_table = base+'.stats'+ending

#load OTU table into pandas DataFrame
df = pd.read_csv(args.otu_table, sep='\t')
df.set_index('OTUId', inplace=True)

ufitslib.log.info("OTU table contains %i OTUs" % len(df.index))

if args.mock_barcode: #if user passes a column name for mock
    #get default mock community value
    if args.mc == "mock3":
        mock = os.path.join(parentdir, 'DB', 'ufits_mock3.fa')
    elif args.mc == "mock2":
        mock = os.path.join(parentdir, 'DB', 'ufits_mock2.fa')
    elif args.mc == "mock1":
        mock = os.path.join(parentdir, 'DB', 'ufits_mock1.fa')
    elif args.mc == "synmock":
        mock = os.path.join(parentdir, 'DB', 'ufits_synmock.fa')
    else:
        mock = os.path.abspath(args.mc)

    #open mock community fasta and count records
    mock_ref_count = ufitslib.countfasta(mock)
    
    #map OTUs to mock community
    mock_out = base + '.mockmap.uc'
    ufitslib.log.info("Mapping OTUs to Mock Community (USEARCH8)")
    if args.mc == 'synmock':
        ufitslib.log.debug("%s -usearch_global %s -strand plus -id 0.80 -db %s -uc %s" % (usearch, args.fasta, mock, mock_out))
        subprocess.call([usearch, '-usearch_global', args.fasta, '-strand', 'plus', '-id', '0.80', '-db', mock, '-uc', mock_out], stdout = FNULL, stderr = FNULL)
    else:
        ufitslib.log.debug("%s -usearch_global %s -strand plus -id 0.97 -db %s -uc %s" % (usearch, args.fasta, mock, mock_out))
        subprocess.call([usearch, '-usearch_global', args.fasta, '-strand', 'plus', '-id', '0.97', '-db', mock, '-uc', mock_out], stdout = FNULL, stderr = FNULL)
    
    #generate dictionary for name change
    annotate_dict = {}
    with open(mock_out, 'rU') as map:
        map_csv = csv.reader(map, delimiter='\t')
        for line in map_csv:
            if line[-1] != "*":
                if float(line[3]) < 97.0:
                    ID = line[-1]+'_spurious_' + line[-2]
                    annotate_dict[line[-2]] = ID
                elif float(line[3]) < 99.0:
                    ID = line[-1]+'_noisy_' + line[-2]
                    annotate_dict[line[-2]] = ID
                elif float(line[3]) < 100.0:
                    ID = line[-1]+'_good_' + line[-2]
                    annotate_dict[line[-2]] = ID                               
                else:
                    ID = line[-1]+'_perfect_'+line[-2]
                    annotate_dict[line[-2]] = ID
else:
    otu_new = args.fasta

#rename OTUs
if args.mock_barcode:
    df.rename(index=annotate_dict, inplace=True)

#sort the table
df2 = df.reindex(index=natsorted(df.index))
if args.col_order == 'naturally':
    ufitslib.log.info("Sorting OTU table naturally")
    df = df2.reindex(columns=natsorted(df2.columns))
else:
    ufitslib.log.info("Sorting OTU table by user defined order (--col_order)")
    col_headers = args.col_order.split(',')
    #check if all names in headers or not
    for i in col_headers:
        if not i in df2.columns.values:
            col_headers.remove(i)
    df = df2.reindex(columns=col_headers)

df.to_csv(sorted_table, sep=delim)

#get sums of columns
fs = df.sum(axis=0)
#fs.to_csv('reads.per.sample.csv')
otus_per_sample_original = df[df > 0].count(axis=0, numeric_only=True)
filtered = pd.DataFrame(df, columns=fs.index)
filt2 = filtered.loc[(filtered != 0).any(1)]
os = filt2.sum(axis=1)
ufitslib.log.info("Removing singleton OTUs (OTUs with only 1 read from all samples)")
fotus = os[os > 2] #valid allele must be found more than 2 times, i.e. no singletons
filt3 = pd.DataFrame(filt2, index=fotus.index)

if args.normalize == 'y':
    #normalize the OTU table
    normal = filt3.truediv(fs)
    normal.to_csv(normal_table_pct, sep=delim)
    #normalize back to read counts, pretend 100,000 reads in each
    norm_round = np.round(normal.multiply(100000), decimals=0)
    norm_round.to_csv(normal_table_nums, sep=delim)
    ufitslib.log.info("Normalizing OTU table to number of reads per sample")
else:
    norm_round = filt3
    
if args.mock_barcode:
    #now calculate the index-bleed in both directions (into the mock and mock into the other samples)
    mock = []
    sample = []
    for i in norm_round.index:
        if not i.startswith('OTU'):
            mock.append(i)
        else:
            sample.append(i)
    mock_df = pd.DataFrame(norm_round, index=mock)
    total = np.sum(np.sum(mock_df,axis=None))
    mock_df.drop(args.mock_barcode, axis=1, inplace=True)
    bleed1 = np.sum(np.sum(mock_df,axis=None))
    bleed1max = bleed1 / float(total)
    sample_df = pd.DataFrame(norm_round, index=sample, columns=[args.mock_barcode])
    bleed2 = np.sum(np.sum(sample_df,axis=None))
    mock_sample = pd.DataFrame(norm_round, columns=[args.mock_barcode])
    bleed2max = bleed2 / float(np.sum(mock_sample.sum(axis=1)))
    
    subtract_num = max(sample_df.max())

    #get max values for bleed
    #can only use into samples measurement if using synmock
    if args.mc == 'synmock':
        if bleed1max > bleed2max:
            bleedfilter = math.ceil(bleed1max*1000)/1000
        else:
            bleedfilter = math.ceil(bleed2max*1000)/1000
        ufitslib.log.info("Index bleed, mock into samples: %f%%.  Index bleed, samples into mock: %f%%." % (bleed1max*100, bleed2max*100))
    else:
        bleedfilter = math.ceil(bleed2max*1000)/1000
        ufitslib.log.info("Index bleed, samples into mock: %f%%." % (bleed2max*100))
        
else:
    bleedfilter = args.index_bleed #this is value needed to filter MiSeq, Ion is likely less, but shouldn't effect the data very much either way.

if args.index_bleed:
    ufitslib.log.info("Overwriting auto detect index-bleed, setting to %f%%" % (args.index_bleed*100))
    bleedfilter = args.index_bleed
else:
    if bleedfilter:
        ufitslib.log.info("Will use value of %f%% for index-bleed OTU filtering." % (bleedfilter*100))
    else:
        bleedfilter = 0 #no filtering if you don't pass -p or -b 
        ufitslib.log.info("No spike-in mock (-b) or index-bleed (-p) specified, thus not running index-bleed filtering") 

#to combat barcode switching, loop through each OTU filtering out if less than bleedfilter threshold
cleaned = []
for row in norm_round.itertuples():
    result = [row[0]]
    total = sum(row[1:])
    sub = total * bleedfilter
    for i in row[1:]:
        if i < sub:
            i = 0
        result.append(i)
    cleaned.append(result)

header = ['OTUId']
for i in norm_round.columns:
    header.append(i)
final = pd.DataFrame(cleaned, columns=header)
final.set_index('OTUId', inplace=True)
if args.subtract != 0:
    if args.subtract != 'auto':
        subtract_num = int(args.subtract)  
    ufitslib.log.info("Subtracting %i from OTU table" % subtract_num)
    sub = final.subtract(subtract_num)
    sub[sub < 0] = 0 #if negative, change to zero
    if not args.keep_mock:
        try:
            sub.drop(args.mock_barcode, axis=1, inplace=True)
        except:
            pass
    sub = sub.loc[~(sub==0).all(axis=1)]
    sub = sub.astype(int)
    sub.to_csv(subtract_table, sep=delim)
    otus_if_sub = sub[sub > 0].count(axis=0, numeric_only=True)
    final = sub.astype(int)
otus_per_sample = final[final > 0].count(axis=0, numeric_only=True)
stats = pd.concat([fs, otus_per_sample_original, otus_per_sample], axis=1)
stats.columns = ['reads per sample', 'original OTUs', 'final OTUs']
print stats.to_string()
stats.to_csv(stats_table, sep=delim)
if not args.keep_mock:
    try:
        final.drop(args.mock_barcode, axis=1, inplace=True)
    except:
        pass
final = final.loc[~(final==0).all(axis=1)]
final = final.astype(int)
final.to_csv(final_table, sep=delim)
final[final > 0] = 1
final.to_csv(final_binary_table, sep=delim)
ufitslib.log.info("Filtering OTU table down to %i OTUs" % (len(final.index)))

#generate final OTU list for taxonomy
ufitslib.log.info("Filtering valid OTUs")
otu_new = base + '.filtered.otus.fa'
with open(otu_new, 'w') as otu_update:
    with open(args.fasta, "rU") as myfasta:
        for rec in SeqIO.parse(myfasta, 'fasta'):
            if args.mock_barcode:
                #map new names of mock
                if rec.id in annotate_dict:
                    newname = annotate_dict.get(rec.id)
                    rec.id = newname
                    rec.description = ''
            if rec.id in final.index:
                SeqIO.write(rec, otu_update, 'fasta')
                
#tell user what output files are
print "-------------------------------------------------------"
print "OTU Table filtering finished"
print "-------------------------------------------------------"
print "OTU Table Stats:   %s" % stats_table
print "Sorted OTU table:  %s" % sorted_table
print "Normalized (pct):  %s" % normal_table_pct
print "Normalized (10k):  %s" % normal_table_nums
if args.subtract != 0:
    print "Subtracted table:  %s" % subtract_table
print "Final filtered:    %s" % final_table
print "Final binary:      %s" % final_binary_table
print "-------------------------------------------------------"

if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits taxonomy -i %s --append_taxonomy %s \n" % (otu_new, final_binary_table)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " ufits taxonomy -i %s --append_taxonomy %s \n" % (otu_new, final_binary_table)
