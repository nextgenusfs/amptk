#!/usr/bin/env python

#This script is a wrapper for -fastq_mergepairs from USEARCH8
from __future__ import print_function
from __future__ import division
from builtins import str
from past.utils import old_div
from builtins import object
import os, sys, argparse, shutil, subprocess, gzip, logging
    
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

class col(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-merge_illumina.py', usage="%(prog)s [options] file_R1.fastq file_R2.fastq",
    description='''Wrapper script for USEARCH8 to merge Illumina PE reads.  Gzip files supported.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-f','--for', dest='fastq_forward_reads', required=True, help='FASTQ R1 file')
parser.add_argument('-r','--rev', dest='fastq_reverse_reads', required=True, help='FASTQ R2 file')
parser.add_argument('-o','--out', dest="out", default='out', help='BaseName for output files')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

def myround(x, base=10):
    return int(base * round(old_div(float(x),base)))

def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    else:
        stdoutformat = logging.Formatter(col.GRN+'%(asctime)s'+col.END+': %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
    log = logging.getLogger("myapp")
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)

#remove logfile if exists
log_name = args.out + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print("-------------------------------------------------------")
#initialize script, log system info and usearch version
log.info("Operating system: %s" % sys.platform)
usearch = args.usearch
try:
    usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
except OSError:
    log.warning("%s not found in your PATH, exiting." % usearch)
    os._exit(1)
log.info("USEARCH version: %s" % usearch_test)


#check extension and decompress if ending in .gz
extension = os.path.splitext(args.fastq_forward_reads)[1]
f1_name = os.path.splitext(args.fastq_forward_reads)[0]
f2_name = os.path.splitext(args.fastq_reverse_reads)[0]
if extension == ".gz":
    log.info("Extracting compressed input files")
    R1_file = gzip.open(args.fastq_forward_reads, 'rb')
    forward = R1_file.read()
    R1_file.close()
    R1_tmp = file(f1_name, 'wb')
    R1_tmp.write(forward)
    R1_tmp.close()
    for_reads = f1_name
    R2_file = gzip.open(args.fastq_reverse_reads, 'rb')
    rev = R2_file.read()
    R2_file.close()
    R2_tmp = file(f2_name, 'wb')
    R2_tmp.write(rev)
    R2_tmp.close()
    rev_reads = f2_name
else:
    for_reads = args.fastq_forward_reads
    rev_reads = args.fastq_reverse_reads

#get read length
fp = open(for_reads)
for i, line in enumerate(fp):
    if i == 1:
        read_length = len(line)
        read_length = myround(read_length)
    elif i > 2:
        break
fp.close()

log.info("Measured Read Length is: %i" % read_length)

#now trim the last bp off of the Illumina data (there for phasing, i.e. 250 bp reads are 251 bp)
pretrim_R1 = args.out + 'pretrim_R1.fq'
pretrim_R2 = args.out + 'pretrim_R2.fq'
log.info("Pre-Processing Reads")
log.debug("%s -fastq_filter %s -fastq_trunclen %s -fastqout %s\n" % (usearch, for_reads, str(read_length), pretrim_R1))
log.debug("%s -fastq_filter %s -fastq_trunclen %s -fastqout %s\n" % (usearch, rev_reads, str(read_length), pretrim_R2))
subprocess.call([usearch, '-fastq_filter', for_reads, '-fastq_trunclen', str(read_length), '-fastqout', pretrim_R1], stdout = FNULL, stderr = FNULL)
subprocess.call([usearch, '-fastq_filter', rev_reads, '-fastq_trunclen', str(read_length), '-fastqout', pretrim_R2], stdout = FNULL, stderr = FNULL)

#next run USEARCH8 mergepe
merge_out = args.out + 'merged.fq'
skip_for = args.out + 'notmerged.R1.fq'
skip_rev = args.out + 'notmerged.R2.fq'
log.info("Merging Overlapping Pairs")
log.debug("%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastqout_notmerged_fwd %s -fastqout_notmerged_rev %s -fastq_truncqual 5 -fastq_allowmergestagger -minhsp 12\n" % (usearch, pretrim_R1, pretrim_R2, merge_out, skip_for, skip_rev))
subprocess.call([usearch, '-fastq_mergepairs', for_reads, '-reverse', rev_reads, '-fastqout', merge_out, '-fastqout_notmerged_fwd', skip_for, '-fastqout_notmerged_rev', skip_rev, '-fastq_truncqual', '5','-fastq_allowmergestagger','-minhsp', '12'], stdout = FNULL, stderr = FNULL)

#now concatenate files for downstream pre-process_illumina.py script
log.info("Concatenating output files")
final_out = args.out + '.fq'
out_file = open(final_out, 'wb')
shutil.copyfileobj(open(merge_out,'rb'), out_file)
shutil.copyfileobj(open(skip_for,'rb'), out_file)
out_file.close()
print("-------------------------------------------------------")
#clean and close up intermediate files
os.remove(merge_out)
os.remove(pretrim_R1)
os.remove(pretrim_R2)
os.remove(skip_for)
os.remove(skip_rev)
