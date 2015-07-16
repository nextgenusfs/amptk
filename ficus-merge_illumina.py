#!/usr/bin/env python

import os
import argparse
import shutil
import subprocess
import gzip
    
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
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

#open log file for usearch8 stderr redirect
log_name = args.out + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)
log_file = open(log_name, 'ab')

#check extension and decompress if ending in .gz
extension = os.path.splitext(args.fastq_forward_reads)[1]
f1_name = os.path.splitext(args.fastq_forward_reads)[0]
f2_name = os.path.splitext(args.fastq_reverse_reads)[0]
if extension == ".gz":
    print "\nExtracting compressed input files"
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

usearch = args.usearch
try:
    subprocess.call([usearch, '--version'], stdout = log_file, stderr = log_file)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)
    
#First run USEARCH8 mergepe
merge_out = args.out + 'merged.fq'
skip_for = args.out + 'notmerged.R1.fq'
skip_rev = args.out + 'notmerged.R2.fq'
print "\nCMD: %s -fastq_mergepairs %s -reverse %s -fastqout %s -fastqout_notmerged_fwd %s -fastqout_notmerged_rev %s\n" % (usearch, for_reads, rev_reads, merge_out, skip_for, skip_rev)
subprocess.call([usearch, '-fastq_mergepairs', for_reads, '-reverse', rev_reads, '-fastqout', merge_out, '-fastqout_notmerged_fwd', skip_for, '-fastqout_notmerged_rev', skip_rev], stdout = log_file, stderr = log_file)

#now join the not merged PE files
join_out = args.out + 'join.fq'
print "CMD: %s -fastq_join %s -reverse %s -fastqout %s\n" % (usearch, skip_for, skip_rev, join_out)
subprocess.call([usearch, '-fastq_join', skip_for, '-reverse', skip_rev, '-fastqout', join_out], stdout = log_file, stderr = log_file)

#now concatenate files for downstream pre-process_illumina.py script
print "Concatenating output files\n"
final_out = args.out + '.fq'
out_file = open(final_out, 'wb')
shutil.copyfileobj(open(merge_out,'rb'), out_file)
shutil.copyfileobj(open(join_out,'rb'), out_file)
out_file.close()

#clean and close up intermediate files
log_file.close()
os.remove(merge_out)
os.remove(skip_for)
os.remove(skip_rev)
os.remove(join_out)
print "------------------------------------------------"
print "Script Finished Successfully!"
print "------------------------------------------------"
print "Merged/Joined FASTQ output:  %s\n" % (final_out)