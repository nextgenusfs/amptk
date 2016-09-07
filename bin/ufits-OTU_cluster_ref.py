#!/usr/bin/env python

#This script runs reference based OTU clustering
#written by Jon Palmer nextgenusfs@gmail.com

import sys, os, argparse, subprocess, inspect, csv, re, logging, shutil, multiprocessing
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib
import pandas as pd

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(script_path)

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

parser=argparse.ArgumentParser(prog='ufits-OTU_cluster_ref.py', usage="%(prog)s [options] -i file.demux.fq\n%(prog)s -h for help menu",
    description='''Script runs UPARSE OTU clustering.
    Requires USEARCH by Robert C. Edgar: http://drive5.com/usearch''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
parser.add_argument('--id', default='97', help="Threshold for alignment")
parser.add_argument('-m','--minsize', default='2', help='Min identical seqs to process')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--uchime_ref', help='Run UCHIME REF [ITS,16S,LSU,COI,custom]')
parser.add_argument('--map_filtered', action='store_true', help='map quality filtered reads back to OTUs')
parser.add_argument('-d','--db', default='USEARCH', help='Reference Database (fasta format)')
parser.add_argument('--utax_db', default='ITS2', help='UTAX Reference Database')
parser.add_argument('--utax_cutoff', default=0.8, type=restricted_float, help='UTAX confidence value threshold.')
parser.add_argument('--utax_level', default='k', choices=['k','p','c','o','f','g','s'], help='UTAX classification level to retain')
parser.add_argument('--mock', default='synmock', help='Spike-in mock community (fasta)')
parser.add_argument('--cleanup', action='store_true', help='Remove Intermediate Files')
parser.add_argument('--closed_ref_only', action='store_true', help='Only run closed reference clustering')
args=parser.parse_args()

def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize

taxonomyLookup = {'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class', 'o': 'Order', 'f': 'Family', 'g': 'Genus', 's': 'Species'}

#remove logfile if exists
log_name = args.out + '.ufits-cluster_ref.log'
if os.path.isfile(log_name):
    os.remove(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
ufitslib.log.info("Operating system: %s, %s" % (sys.platform, ufitslib.get_version()))
usearch = args.usearch
try:
    usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
except OSError:
    ufitslib.log.warning("%s not found in your PATH, exiting." % usearch)
    os._exit(1)
ufitslib.log.info("USEARCH version: %s" % usearch_test)

#check if vsearch version > 1.9.1 is installed
vsearch_check = ufitslib.which('vsearch')
if vsearch_check:
    vsearch = ufitslib.checkvsearch()
    vsearch_version = ufitslib.get_vsearch_version()
    if vsearch:
        ufitslib.log.info("vsearch v%s detected, will use for filtering" % vsearch_version)
    else:
        ufitslib.log.info("vsearch v%s detected, need version at least v1.9.1, using Python for filtering")
else:
    vsearch = False
    ufitslib.log.info("vsearch not installed, using Python for filtering")

#make tmp folder
tmp = args.out + '_tmp'
if not os.path.exists(tmp):
    os.makedirs(tmp)

#get utax_database
installed_utax = []
installed_db = []
for file in os.listdir(os.path.join(parentdir, 'DB')):
    if file.endswith('.udb.txt'):
        filename = os.path.join(parentdir, 'DB', file)
        with open(filename) as input:
            for line in input:
                if line.startswith('utax'):
                    installed_utax.append(file.replace('.udb.txt',''))  
    elif file.endswith('.extracted.fa'):
        installed_db.append(file)             
if args.utax_db in installed_utax:
    ufitslib.log.info("Using %s UTAX database" % args.utax_db)
    utaxDB = os.path.join(parentdir, 'DB', args.utax_db+'.udb')
else:
    ufitslib.log.info("Custom UTAX database entered, you can install one with `ufits database` command.")
    utaxDB = os.abspath(args.utax_db)

#format reference database and look for duplicate annotation, error if found
ufitslib.log.info("Checking Reference Database")
installed_prefix = []
for i in installed_db:
    installed_prefix.append(i.split('.')[0])
if args.db in installed_prefix:
    match = installed_prefix.index(args.db)
    DB = os.path.join(parentdir, 'DB', installed_db[match])
else:
    DB = os.path.abspath(args.db)
refDB = os.path.join(tmp, 'reference_DB.fa')
if args.mock:
    if args.mock == 'synmock':
        mock = os.path.join(parentdir, 'DB', 'ufits_synmock.fa')
    else:
        mock = os.path.abspath(args.mock)
seen = []
with open(refDB, 'w') as output:
    if args.mock:
        with open(mock) as input1:
            for rec in SeqIO.parse(input1, 'fasta'):
                if not rec.id in seen:
                    SeqIO.write(rec, output, 'fasta')
                else:
                    ufitslib.log.error("Duplicate ID's in Ref DB: %s, exiting" % rec.id)
                    os._exit(1)
    with open(DB) as input2:
        for rec in SeqIO.parse(input2, 'fasta'):
            if not rec.id in seen:
                SeqIO.write(rec, output, 'fasta')
            else:
                ufitslib.log.error("Duplicate ID's in Ref DB: %s, exiting" % rec.id)
                os._exit(1)

#Count FASTQ records
ufitslib.log.info("Loading FASTQ Records")
orig_total = ufitslib.countfastq(args.FASTQ)
size = checkfastqsize(args.FASTQ)
readablesize = ufitslib.convertSize(size)
ufitslib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

#Expected Errors filtering step and convert to fasta
filter_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.filter.fq')
filter_fasta = os.path.join(tmp, args.out + '.EE' + args.maxee + '.filter.fa')
orig_fasta = os.path.join(tmp, args.out+'.orig.fa')
ufitslib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
if vsearch:
    subprocess.call(['vsearch', '--fastq_filter', args.FASTQ, '--fastq_maxee', str(args.maxee), '--fastqout', filter_out, '--fastaout', filter_fasta, '--fastq_qmax', '45'], stdout = FNULL, stderr = FNULL)
    subprocess.call(['vsearch', '--fastq_filter', args.FASTQ, '--fastaout', orig_fasta, '--fastq_qmax', '45'], stdout = FNULL, stderr = FNULL)
else:
    with open(filter_out, 'w') as output:
        SeqIO.write(ufitslib.MaxEEFilter(args.FASTQ, args.maxee), output, 'fastq')
    SeqIO.convert(args.FASTQ, 'fastq', orig_fasta, 'fasta')
    SeqIO.convert(filter_out, 'fastq', filter_fasta, 'fasta')
qtrimtotal = ufitslib.countfastq(filter_out)
ufitslib.log.info('{0:,}'.format(qtrimtotal) + ' reads passed')

#now run full length dereplication
derep_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.derep.fa')
ufitslib.log.info("De-replication (remove duplicate reads)")
if vsearch:
    subprocess.call(['vsearch', '--derep_fulllength', filter_fasta, '--sizeout', '--output', derep_out], stdout = FNULL, stderr = FNULL)
else:
    ufitslib.dereplicate(filter_out, derep_out)
total = ufitslib.countfasta(derep_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#now run usearch 8 sort by size
sort_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.sort.fa')
ufitslib.log.info("Sorting reads by size: removing reads seen less than %s times" % args.minsize)
if vsearch:
    ufitslib.log.debug("vsearch --sortbysize %s --minsize %s --output %s" % (derep_out, args.minsize, sort_out))
    subprocess.call(['vsearch', '--sortbysize', derep_out, '--minsize', args.minsize, '--output', sort_out], stdout = FNULL, stderr = FNULL)
else:
    ufitslib.log.debug("%s -sortbysize %s -minsize %s -fastaout %s" % (usearch, unoise_out, args.minsize, sort_out))
    subprocess.call([usearch, '-sortbysize', unoise_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = FNULL, stderr = FNULL)
total = ufitslib.countfasta(sort_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#chimera detection
#first run through de novo chimera detection, using cluster_otus as that is apparently recommended over uchime_denovo
ufitslib.log.info("De novo chimera detection (USEARCH)")
chimera_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.chimera_check.fa')
ufitslib.log.debug("%s -cluster_otus %s -sizein -sizeout -relabel OTU_ -otu_radius_pct 0 -otus %s" % (usearch, sort_out, chimera_out))
subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-sizeout', '-relabel', 'chimera_', '-otu_radius_pct', '0', '-otus', chimera_out], stdout = FNULL, stderr = FNULL)
total = ufitslib.countfasta(chimera_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')
  
#optional UCHIME Ref
if not args.uchime_ref:
    uchime_out = chimera_out
else:
    uchime_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.uchime.otus.fa')
    #R. Edgar now says using largest DB is better for UCHIME, so use the one distributed with taxonomy
    if args.uchime_ref in ['ITS', '16S', 'LSU', 'COI']: #test if it is one that is setup, otherwise default to full path
        uchime_db = os.path.join(parentdir, 'DB', args.uchime_ref+'.extracted.fa')
        if not os.path.isfile(uchime_db):
            ufitslib.error("Database not properly configured, run `ufits install` to setup DB, skipping chimera filtering")
            uchime_out = chimera_out
    else:
        uchime_db = os.path.abspath(args.uchime_ref)
    #now run chimera filtering if all checks out
    if not os.path.isfile(uchime_out):
        if vsearch:
            ufitslib.log.info("Chimera Filtering (VSEARCH)")
            ufitslib.log.debug("vsearch --uchime_ref %s --db %s --nonchimeras %s --mindiv 1" % (chimera_out, uchime_db, uchime_out))
            subprocess.call(['vsearch', '--mindiv', '1.0', '--uchime_ref', otu_clean, '--db', uchime_db, '--nonchimeras', uchime_out], stdout = FNULL, stderr = FNULL)
        else:
            ufitslib.log.info("Chimera Filtering (UCHIME)")
            ufitslib.log.debug("%s -uchime_ref %s -strand plus -db %s -nonchimeras %s -mindiv 1" % (usearch, chimera_out, uchime_db, uchime_out))
            subprocess.call([usearch, '-uchime_ref', chimera_out, '-strand', 'plus', '-db', uchime_db, '-nonchimeras', uchime_out, '-mindiv', '1.0'], stdout = FNULL, stderr = FNULL)
        total = ufitslib.countfasta(uchime_out)
        ufitslib.log.info('{0:,}'.format(total) + ' OTUs passed')
    
    
#now run usearch_global versus reference database
align_out = os.path.join(tmp, args.out + '.align.uc')
pident = int(args.id) * 0.01
ufitslib.log.info("Reference Clustering using Global Alignment, %s%% identity" % args.id)
if vsearch:
    subprocess.call(['vsearch', '--usearch_global', uchime_out, '--db', refDB, '--id', str(pident), '--output_no_hits', '--top_hits_only', '--notrunclabels', '--uc', align_out], stdout=FNULL, stderr=FNULL)
else:
     subprocess.call([usearch, '-usearch_global', uchime_out, '-db', refDB, '-id', str(pident), '-strand', 'plus', '--top_hit_only', '--uc', align_out, '-notrunclabels'], stdout=FNULL, stderr=FNULL)

#parse results
ref_results = {}
nohits = []
with open(align_out, 'rU') as alignment:
    for line in alignment:
        line = line.replace('\n', '')
        col = line.split('\t')
        counts = col[8].split(';')
        counts = int(counts[1].replace('size=', ''))
        if col[3] == '*':
            nohits.append(col[8])
            continue
        if float(col[3]) >= float(args.id):
            if not col[8] in ref_results:
                ref_results[col[8]] = (col[9], col[3], counts)
            else:
                print "Error: %s duplicated ID" % col[8]
        else:
            nohits.append(col[8])

#summarize results from first ref clustering
num_refcluster = len(ref_results)
seqs_refcluster = 0
for k,v in ref_results.items():
    seqs_refcluster += v[2]
ufitslib.log.info("%i OTUs classified " % num_refcluster + "({0:.0f}%".format(seqs_refcluster/float(qtrimtotal)* 100)+ " of reads)")

#get ref clustered hits to file with taxonomy
ref_clustered = os.path.join(tmp, args.out+'.ref_clustered.fa')
with open(ref_clustered, 'w') as refoutput:
    with open(uchime_out, 'rU') as input:
        otu_counter = 1
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id in ref_results:
                res = ref_results.get(rec.id)
                pident = res[1]
                tax = res[0]
                newID = 'OTU_'+str(otu_counter)+';pident='+pident+';'+tax
                rec.id = newID
                rec.name = ''
                rec.description = ''
                SeqIO.write(rec, refoutput, 'fasta')
                otu_counter += 1

if not args.closed_ref_only:
    #get nohits file to run clustering
    utax_ref = os.path.join(tmp, args.out + '.EE' + args.maxee + '.utax_ref.fa')
    with open(utax_ref, 'w') as output:
        with open(uchime_out, 'rU') as input:
            for rec in SeqIO.parse(input, 'fasta'):
                if rec.id in nohits:
                    SeqIO.write(rec, output, 'fasta')

    #input needs to be sorted, so 
    ref_sort = os.path.join(tmp, args.out+'.utax_ref.sorted.fa')
    if vsearch:
        ufitslib.log.debug("vsearch --sortbysize %s --minsize %s --output %s" % (utax_ref, args.minsize, ref_sort))
        subprocess.call(['vsearch', '--sortbysize', utax_ref, '--minsize', args.minsize, '--output', ref_sort], stdout = FNULL, stderr = FNULL)
    else:
        ufitslib.log.debug("%s -sortbysize %s -minsize %s -fastaout %s" % (usearch, utax_ref, args.minsize, ref_sort))
        subprocess.call([usearch, '-sortbysize', utax_ref, '-minsize', args.minsize, '-fastaout', ref_sort], stdout = FNULL, stderr = FNULL)
           
    #now run clustering algorithm on those not found in reference database
    radius = str(100 - int(args.pct_otu))
    otu_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.otus.fa')
    ufitslib.log.info("De novo Clustering remaining sequences (UPARSE)")
    ufitslib.log.debug("%s -cluster_otus %s -sizein -sizeout -relabel OTU_ -otu_radius_pct %s -otus %s" % (usearch, ref_sort, radius, otu_out))
    subprocess.call([usearch, '-cluster_otus', ref_sort, '-sizein', '-sizeout', '-relabel', 'OTU_', '-otu_radius_pct', radius, '-otus', otu_out], stdout=FNULL, stderr=FNULL)
    total = ufitslib.countfasta(otu_out)
    ufitslib.log.info('{0:,}'.format(total) + ' de novo OTUs')

    #try utax reference clustering
    ufitslib.log.info("Reference Clustering de novo OTUs using UTAX")
    subprocess.call([usearch, '-cluster_otus_utax', otu_out, '-db', utaxDB, '-utax_cutoff', str(args.utax_cutoff), '-utax_level', 's', '-strand', 'plus', '-utaxout', os.path.join(tmp, args.out+'.utax.out')], stdout=FNULL, stderr=FNULL)
    #setup tax filtering
    tax_values = ['k','p','c','o','f','g','s']
    filter_index = tax_values.index(args.utax_level)
    filt_tax_values = [s + ':' for s in tax_values[filter_index:]]
    #get results from utax
    with open(ref_clustered, 'a') as output:
        seqDict = SeqIO.index(otu_out, 'fasta')
        utaxresults = []
        with open(os.path.join(tmp, args.out+'.utax.out'), 'ru') as utax:
            for line in utax:
                line = line.replace('\n', '')
                col = line.split('\t')
                ID = col[0]
                tax = col[2]
                if any(x in tax for x in filt_tax_values):
                    record = seqDict[ID]
                    record.id = 'OTU_'+str(otu_counter)+';UTAX;tax='+tax
                    record.name = ''
                    record.description = ''
                    SeqIO.write(record, output, 'fasta')
                    otu_counter += 1
    total = ufitslib.countfasta(ref_clustered) - num_refcluster
    ufitslib.log.info('{0:,}'.format(total) + ' classified to %s' % taxonomyLookup.get(args.utax_level))

#clean up padded N's
ufitslib.log.info("Cleaning up padding from OTUs")
otu_clean = os.path.join(tmp, args.out + '.clean.otus.fa')
with open(otu_clean, 'w') as output:
    with open(ref_clustered, 'rU') as input:
        for line in input:
            if line.startswith (">"):
                output.write(line)
            else:
                line = re.sub('[^GATC]', "", line)
                line = line + '\n'
                if line != '\n':
                    output.write(line)             
total = ufitslib.countfasta(otu_clean)
ufitslib.log.info('{0:,}'.format(total) + ' total OTUs')
       
#now map reads back to OTUs
uc_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.mapping.uc')
if args.map_filtered:
    reads = filter_fasta
else:
    reads = orig_fasta

ufitslib.log.info("Mapping Reads to OTUs")
if vsearch:
    subprocess.call(['vsearch', '-usearch_global', reads, '-strand', 'plus', '-id', '0.97', '-db', otu_clean, '-uc', uc_out, '--notrunclabels'], stdout = FNULL, stderr = FNULL)
else:
    ufitslib.log.debug("%s -usearch_global %s -strand plus -id 0.97 -db %s -uc %s" % (usearch, reads, otu_clean, uc_out))
    subprocess.call([usearch, '-usearch_global', reads, '-strand', 'plus', '-id', '0.97', '-db', otu_clean, '-uc', uc_out, '-notrunclabels'], stdout = FNULL, stderr = FNULL)

#count reads mapped
if vsearch:
    total = ufitslib.line_count(uc_out)
else:
    total = ufitslib.line_count2(uc_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

#Build OTU table
otu_table = os.path.join(tmp, args.out + '.EE' + args.maxee + '.otu_table.txt')
uc2tab = os.path.join(parentdir, 'lib', 'uc2otutable.py')
ufitslib.log.info("Creating OTU Table")
ufitslib.log.debug("%s %s %s" % (uc2tab, uc_out, otu_table))
subprocess.call([sys.executable, uc2tab, uc_out, otu_table], stdout = FNULL, stderr = FNULL)

#Move files around, delete tmp if argument passed.
currentdir = os.getcwd()
final_otu = os.path.join(currentdir, args.out + '.cluster.otus.fa')
shutil.copyfile(otu_clean, final_otu)
final_otu_table = os.path.join(currentdir, args.out + '.otu_table.txt')
#shutil.copyfile(otu_table, final_otu_table)

#split the taxonomy from otu name, making new column at the end
with open(final_otu_table, 'w') as output:
    with open(otu_table, 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            cols = line.split('\t')
            if line.startswith('OTUId'):
                cols.append('Taxonomy')
            else:
                otuname = cols[0].split(';', 1)[0]
                tax = cols[0].split(';', 1)[1]
                cols.append(tax)
                cols[0] = cols[0].replace(';'+tax, '')
            output.write('%s\n' % '\t'.join(cols))

if args.cleanup:
    shutil.rmtree(tmp)

#Print location of files to STDOUT
print "-------------------------------------------------------"
print "OTU Clustering Script has Finished Successfully"
print "-------------------------------------------------------"
if not args.cleanup:
    print "Tmp Folder of files: %s" % tmp
print "Clustered OTUs: %s" % final_otu
print "OTU Table: %s" % final_otu_table
print "-------------------------------------------------------"

otu_print = final_otu.split('/')[-1]
tab_print = final_otu_table.split('/')[-1]
if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END + " ufits filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print)
