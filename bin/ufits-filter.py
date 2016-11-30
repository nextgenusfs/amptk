#!/usr/bin/env python

#This script filters results from ufits-OTU_cluster.py
#written by Jon Palmer palmer.jona at gmail dot com

import sys, os, argparse, inspect, subprocess, csv, sys, math
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
parser.add_argument('-p','--index_bleed',  help='Index Bleed filter. Default: auto')
parser.add_argument('-s','--subtract', default=0, help='Threshold to subtract')
parser.add_argument('-n','--normalize', default='y', choices=['y','n'], help='Normalize OTU table prior to filtering')
parser.add_argument('--mc', help='Multi-FASTA mock community')
parser.add_argument('-d','--delimiter', default='tsv', choices=['csv','tsv'], help='Delimiter')
parser.add_argument('--col_order', dest="col_order", default="naturally", help='Provide comma separated list')
parser.add_argument('--keep_mock', action='store_true', help='Keep mock sample in OTU table (Default: False)')
parser.add_argument('--show_stats', action='store_true', help='Show stats datatable STDOUT')
parser.add_argument('-o','--out', help='Base output name')
parser.add_argument('--min_reads_otu', default=2, type=int, help='Minimum number of reads per OTU for experiment')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH8 EXE')
parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
args=parser.parse_args()

if not args.out:
    #get base name of files
    base = args.otu_table.split(".otu_table")
    base = base[0]
else:
    base = args.out

#remove logfile if exists
log_name = base + '.ufits-filter.log'
ufitslib.removefile(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
ufitslib.SystemInfo()
#Do a version check
usearch = args.usearch
ufitslib.versionDependencyChecks(usearch)

#check if otu_table is empty
ufitslib.log.info("Loading OTU table: %s" % args.otu_table)
check = os.stat(args.otu_table).st_size
if check == 0:
    ufitslib.log.error("Input OTU table is empty")
    sys.exit(1)
#get the OTU header info (depending on how OTU table was constructed, this might be different, so find it as you need for indexing)
with open(args.otu_table, 'rU') as f:
    first_line = f.readline()
    OTUhead = first_line.split('\t')[0]

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
df.set_index(OTUhead, inplace=True)

ufitslib.log.info("OTU table contains %i OTUs" % len(df.index))

#setup output files/variables
mock_out = base + '.mockmap.uc'
mock_sort = base + '.mockmap.sort.uc'

if args.mock_barcode: #if user passes a column name for mock
    #check if mock barcode is valid
    validBCs = df.columns.values.tolist()
    if not args.mock_barcode in validBCs:
        ufitslib.log.error("%s not a valid barcode." % args.mock_barcode)
        ufitslib.log.error("Valid barcodes: %s" % (' '.join(validBCs)))
        sys.exit(1)
    #make sure there is a --mc passed here otherwise throw error
    if not args.mc:
        ufitslib.log.error("If using the -b,--barcode option you must specify a fasta file of mock community via the --mc option")
        sys.exit(1)
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
    ufitslib.log.info("Mapping OTUs to Mock Community (USEARCH)")
    cmd = [usearch, '-usearch_global', mock, '-strand', 'plus', '-id', '0.95', '-db', args.fasta, '-uc', mock_out, '-maxaccepts', '3']
    ufitslib.runSubprocess(cmd, ufitslib.log)
    #sort the output to avoid problems
    with open(mock_sort, 'w') as output:
        subprocess.call(['sort', '-k4,4nr', mock_out], stdout = output)

    #generate dictionary for name change
    found_dict = {}
    missing = []
    chimeras = []
    seen = []
    with open(mock_sort, 'rU') as map:
        map_csv = csv.reader(map, delimiter='\t')
        for line in map_csv:
            if line[-1] != "*":
                if not line[-2] in found_dict:
                    found_dict[line[-2]] = (line[-1], float(line[3]))
                    seen.append(line[-1])
                else:
                    oldpident = found_dict.get(line[-2])[1]
                    oldid = found_dict.get(line[-2])[0]
                    if float(line[3]) > oldpident:
                        found_dict[line[-2]] = (line[-1], float(line[3]))
                        if not oldid in seen:
                            chimeras.append(oldid)
                    else:
                        if not line[-1] in seen:
                            chimeras.append(line[-1])        
            else:
                missing.append(line[-2].split(' ')[0])

    if missing:
        ufitslib.log.info("Mock members not found: %s" % ', '.join(missing))
    #make name change dict
    annotate_dict = {}
    for k,v in found_dict.items():
        ID = v[0].replace('_chimera', '')
        newID = k+'_pident='+str(v[1])+'_'+v[0]
        annotate_dict[ID] = newID
    for i in chimeras:
        annotate_dict[i] = i+'_suspect_mock_chimera'
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
tos = filt2.sum(axis=1)
ufitslib.log.info("Removing OTUs according to --min_reads_otu: (OTUs with less than %i reads from all samples)" % args.min_reads_otu)
fotus = tos[tos >= args.min_reads_otu] #valid allele must be found atleast from than 2 times, i.e. no singletons
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
    #get names from mapping
    for k,v in annotate_dict.items():
        if not '_suspect_mock_chimera' in v:
            mock.append(v)
    for i in norm_round.index:
        if not i in mock:
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
    args.index_bleed = float(args.index_bleed)
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
    total = max(row[1:]) #get max OTU count from table to calculate index bleed from.
    sub = total * bleedfilter
    for i in row[1:]:
        if i < sub:
            i = 0
        result.append(i)
    cleaned.append(result)

header = [OTUhead]
for i in norm_round.columns:
    header.append(i)
final = pd.DataFrame(cleaned, columns=header)
final.set_index(OTUhead, inplace=True)
if args.subtract != 'auto':
    subtract_num = int(args.subtract)
else:
    try:
        subtract_num = int(subtract_num)
        ufitslib.log.info("Auto subtract filter set to %i" % subtract_num)
    except NameError:
        subtract_num = 0
        ufitslib.log.info("Error: to use 'auto' subtract feature, provide a sample name to -b,--mock_barcode.")
if subtract_num != 0:
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
stats.fillna(0, inplace=True)
stats = stats.astype(int)
if args.show_stats:
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
if not args.debug:
    for i in [normal_table_pct, normal_table_nums, subtract_table, mock_out, mock_sort]:
        ufitslib.removefile(i)
else:
    
    print "Normalized (pct):  %s" % normal_table_pct
    print "Normalized (10k):  %s" % normal_table_nums
    if args.subtract != 0:
        print "Subtracted table:  %s" % subtract_table
print "Final filtered:    %s" % final_table
print "Final binary:      %s" % final_binary_table
print "Filtered OTUs:     %s" % otu_new
print "-------------------------------------------------------"

if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits taxonomy -f %s -i %s -m mapping_file.txt -d ITS2\n" % (otu_new, final_binary_table)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " ufits taxonomy -f %s -i %s -m mapping_file.txt -d ITS2\n" % (otu_new, final_binary_table)
