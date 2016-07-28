#!/usr/bin/env python

#This script runs USEARCH OTU clustering
#written by Jon Palmer nextgenusfs@gmail.com

import sys, os, argparse, subprocess, inspect, csv, re, logging, shutil, multiprocessing
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.ufitslib as ufitslib

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

parser=argparse.ArgumentParser(prog='ufits-OTU_cluster.py', usage="%(prog)s [options] -i file.demux.fq\n%(prog)s -h for help menu",
    description='''Script runs UPARSE OTU clustering.
    Requires USEARCH by Robert C. Edgar: http://drive5.com/usearch''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
parser.add_argument('-o','--out', default='out', help='Base output name')
parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
parser.add_argument('-m','--minsize', default='2', help='Min size to keep for clustering')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--uchime_ref', help='Run UCHIME REF [ITS1,ITS2,Full,16S]')
parser.add_argument('--map_filtered', action='store_true', help='map quality filtered reads back to OTUs')
parser.add_argument('--unoise', action='store_true', help='Run De-noising (UNOISE)')
parser.add_argument('--size_annotations', action='store_true', help='Append size annotations')
parser.add_argument('--cleanup', action='store_true', help='Remove Intermediate Files')
args=parser.parse_args()

def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize

#remove logfile if exists
log_name = args.out + '.ufits-cluster.log'
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
total = ufitslib.countfastq(filter_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#now run full length dereplication
derep_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.derep.fa')
ufitslib.log.info("De-replication (remove duplicate reads)")
if vsearch:
    subprocess.call(['vsearch', '--derep_fulllength', filter_fasta, '--sizeout', '--output', derep_out], stdout = FNULL, stderr = FNULL)
else:
    ufitslib.dereplicate(filter_out, derep_out)
total = ufitslib.countfasta(derep_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads passed')

#optional run UNOISE
if args.unoise:
    unoise_out = unoise_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.denoised.fa')
    ufitslib.log.info("Denoising Data with UNOISE")
    ufitslib.log.debug("%s -cluster_fast %s -centroids %s -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size" % (usearch, derep_out, unoise_out))
    subprocess.call([usearch, '-cluster_fast', derep_out, '-centroids', unoise_out, '-id', '0.9', '-maxdiffs', '5', '-abskew', '10', '-sizein', '-sizeout', '-sort', 'size'], stdout = FNULL, stderr = FNULL)
    total = ufitslib.countfasta(unoise_out)
    ufitslib.log.info('{0:,}'.format(total) + ' reads passed')
else:
    unoise_out = derep_out

#now run usearch 8 sort by size
sort_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.sort.fa')
ufitslib.log.info("Sorting reads by size")
ufitslib.log.debug("%s -sortbysize %s -minsize %s -fastaout %s" % (usearch, unoise_out, args.minsize, sort_out))
subprocess.call([usearch, '-sortbysize', unoise_out, '-minsize', args.minsize, '-fastaout', sort_out], stdout = FNULL, stderr = FNULL)


#now run clustering algorithm
radius = str(100 - int(args.pct_otu))
otu_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.otus.fa')
ufitslib.log.info("Clustering OTUs (UPARSE)")
if args.size_annotations:
    ufitslib.log.debug("%s -cluster_otus %s -sizein -sizeout -relabel OTU_ -otu_radius_pct %s -otus %s" % (usearch, sort_out, radius, otu_out))
    subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-sizeout', '-relabel', 'OTU_', '-otu_radius_pct', radius, '-otus', otu_out], stdout = FNULL, stderr = FNULL)
else:
    ufitslib.log.debug("%s -cluster_otus %s -relabel OTU_ -otu_radius_pct %s -otus %s" % (usearch, sort_out, radius, otu_out))
    subprocess.call([usearch, '-cluster_otus', sort_out, '-sizein', '-relabel', 'OTU_', '-otu_radius_pct', radius, '-otus', otu_out], stdout = FNULL, stderr = FNULL)
total = ufitslib.countfasta(otu_out)
ufitslib.log.info('{0:,}'.format(total) + ' OTUs')

#clean up padded N's
ufitslib.log.info("Cleaning up padding from OTUs")
otu_clean = os.path.join(tmp, args.out + '.EE' + args.maxee + '.clean.otus.fa')
with open(otu_clean, 'w') as output:
    with open(otu_out, 'rU') as input:
        for line in input:
            if line.startswith (">"):
                output.write(line)
            else:
                line = re.sub('[^GATC]', "", line)
                line = line + '\n'
                if line != '\n':
                    output.write(line)

#optional UCHIME Ref
if not args.uchime_ref:
    uchime_out = otu_clean
else:
    uchime_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.uchime.otus.fa')
    #You might need to update these in the future, but leaving data and version in name so it is obvious where they came from
    for file in os.listdir(os.path.join(parentdir, 'DB')):
        if file.startswith('uchime_sh'):
            if file.endswith('ITS1.fasta'):
                uchime_ITS1 = os.path.join(parentdir, 'DB', file)
            if file.endswith('ITS2.fasta'):
                uchime_ITS2 = os.path.join(parentdir, 'DB', file)
        if file.startswith('uchime_reference_'):
            uchime_FULL = os.path.join(parentdir, 'DB', file)
    if args.uchime_ref == "ITS1":
        uchime_db = uchime_ITS1
    elif args.uchime_ref == "ITS2":
        uchime_db = uchime_ITS2
    elif args.uchime_ref == "Full":
        uchime_db = uchime_FULL
    elif args.uchime_ref == "16S":
        uchime_db = os.path.join(parentdir, 'DB', 'rdp_gold.fa')
    else:
        uchime_db = os.path.abspath(args.uchime_ref)
    if vsearch:
        ufitslib.log.info("Chimera Filtering (VSEARCH)")
        ufitslib.log.debug("vsearch --uchime_ref %s --db %s --nonchimeras %s --mindiv 1" % (otu_clean, uchime_db, uchime_out))
        subprocess.call(['vsearch --mindiv 1.0 --uchime_ref', otu_clean, '--db', uchime_db, '--nonchimeras', uchime_out])#, stdout = FNULL, stderr = FNULL)
    else:
        ufitslib.log.info("Chimera Filtering (UCHIME): make sure DB has rev_comp'd sequences")
        ufitslib.log.debug("%s -uchime_ref %s -strand plus -db %s -nonchimeras %s -mindiv 1" % (usearch, otu_clean, uchime_db, uchime_out))
        subprocess.call([usearch, '-uchime_ref', otu_clean, '-strand', 'plus', '-db', uchime_db, '-nonchimeras', uchime_out, '-mindiv 1.0'], stdout = FNULL, stderr = FNULL)
    total = ufitslib.countfasta(uchime_out)
    ufitslib.log.info('{0:,}'.format(total) + ' OTUs passed')


#now map reads back to OTUs
uc_out = os.path.join(tmp, args.out + '.EE' + args.maxee + '.mapping.uc')
if args.map_filtered:
    reads = filter_fasta
else:
    reads = orig_fasta

ufitslib.log.info("Mapping Reads to OTUs")
if vsearch:
    subprocess.call(['vsearch', '-usearch_global', reads, '-strand', 'plus', '-id', '0.97', '-db', uchime_out, '-uc', uc_out], stdout = FNULL, stderr = FNULL)
else:
    ufitslib.log.debug("%s -usearch_global %s -strand plus -id 0.97 -db %s -uc %s" % (usearch, reads, uchime_out, uc_out))
    subprocess.call([usearch, '-usearch_global', reads, '-strand', 'plus', '-id', '0.97', '-db', uchime_out, '-uc', uc_out], stdout = FNULL, stderr = FNULL)

#count reads mapped
if vsearch:
    total = ufitslib.line_count(uc_out)
else:
    total = ufitslib.line_count2(uc_out)
ufitslib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs')

#Build OTU table
otu_table = os.path.join(tmp, args.out + '.EE' + args.maxee + '.otu_table.txt')
uc2tab = os.path.join(parentdir, 'lib', 'uc2otutable.py')
ufitslib.log.info("Creating OTU Table")
ufitslib.log.debug("%s %s %s" % (uc2tab, uc_out, otu_table))
subprocess.call([sys.executable, uc2tab, uc_out, otu_table], stdout = FNULL, stderr = FNULL)

#Move files around, delete tmp if argument passed.
currentdir = os.getcwd()
final_otu = os.path.join(currentdir, args.out + '.cluster.otus.fa')
shutil.copyfile(uchime_out, final_otu)
final_otu_table = os.path.join(currentdir, args.out + '.otu_table.txt')
shutil.copyfile(otu_table, final_otu_table)
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
