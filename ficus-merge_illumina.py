#!/usr/bin/env python

import os
import argparse
import shutil
import subprocess
import gzip
import os.path

class bcolors:
    GREEN = '\033[92m'
    BLUE = '\033[36m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

parser=argparse.ArgumentParser(prog='ficus-merge_illumina.py', usage="%(prog)s [options] file_R1.fastq file_R2.fastq",
    description='''Wrapper script for USEARCH8 to merge Illumina PE reads.  Gzip files supported.''',
    epilog="""Written by Jon Palmer (2015) palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('fastq_forward_reads', help='FASTQ R1 file')
parser.add_argument('fastq_reverse_reads', help='FASTQ R2 file')
parser.add_argument('-o','--out_prefix', dest="out", default='out', help='BaseName for output files')
args=parser.parse_args()

extension = os.path.splitext(args.fastq_forward_reads)[1]
f1_name = os.path.splitext(args.fastq_forward_reads)[0]
f2_name = os.path.splitext(args.fastq_reverse_reads)[0]
if extension == ".gz":
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

usearch = "usearch8"
try:
    print "------------------------------------------------"
    print bcolors.BLUE + "Looking for USEARCH8 in your PATH:" + bcolors.ENDC
    subprocess.call([usearch, '--version'])
    print "------------------------------------------------"
except OSError:
    print + bcolors.FAIL + "%s not found in your PATH" % usearch + bcolors.ENDC
    usearch = raw_input("Enter full path to USEARCH or type exit to quit: ")
    
if usearch == "exit":
    os._exit(1)
    
#First run USEARCH8 mergepe
print bcolors.BLUE + "Running USEARCH fastq_mergepairs" + bcolors.ENDC
print "------------------------------------------------"
merge_out = args.out + 'merged.fq'
skip_for = args.out + 'notmerged.R1.fq'
skip_rev = args.out + 'notmerged.R2.fq'
os.system('%s %s %s %s %s %s %s %s %s %s %s' % (usearch, '-fastq_mergepairs', for_reads, '-reverse', rev_reads, '-fastqout', merge_out, '-fastqout_notmerged_fwd', skip_for, '-fastqout_notmerged_rev', skip_rev))

#now join the not merged PE files
print "------------------------------------------------"
print bcolors.BLUE + "Running USEARCH fastq_join" + bcolors.ENDC
print "------------------------------------------------"
join_out = args.out + 'join.fq'
os.system('%s %s %s %s %s %s %s' % (usearch, '-fastq_join', skip_for, '-reverse', skip_rev, '-fastqout', join_out))

#now concatenate files for downstream pre-process_illumina.py script
final_out = args.out + '.fq'
out_file = open(final_out, 'wb')
shutil.copyfileobj(open(merge_out,'rb'), out_file)
shutil.copyfileobj(open(join_out,'rb'), out_file)
out_file.close()

#clean up intermediate files
os.remove(merge_out)
os.remove(skip_for)
os.remove(skip_rev)
os.remove(join_out)
print "------------------------------------------------"
print bcolors.GREEN + "Script finished successfully!" + bcolors.ENDC
print "------------------------------------------------"