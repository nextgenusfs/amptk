#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import inspect
import argparse
import shutil
import logging
import subprocess
import multiprocessing
import glob
import itertools
import re
import edlib
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from natsort import natsorted
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
class col(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def processReadsPE(input, args=False):
    base = os.path.basename(input)
    forward_reads = os.path.join(tmpdir, base+'_R1.fq')
    reverse_reads = os.path.join(tmpdir, base+'_R2.fq')
    index_reads = os.path.join(tmpdir, base+'_R3.fq')
    trim_forward = os.path.join(tmpdir, base+'_R1.trimmed.fq')
    trim_reverse = os.path.join(tmpdir, base+'_R2.trimmed.fq')
    merged_reads = os.path.join(tmpdir, base+'.merged.fq')
    DemuxOut = os.path.join(tmpdir, base+'.demux.fq')
    Total, BCFound, ForPrimerCount, RevPrimerCount = amptklib.DemuxIllumina(forward_reads, reverse_reads, index_reads, Barcodes, args.barcode_mismatch, FwdPrimer, RevPrimer, args.primer_mismatch, trim_forward, trim_reverse)
    amptklib.MergeReadsSimple(trim_forward, trim_reverse, '.', merged_reads, args.min_len, usearch, args.rescue_forward, args.merge_method)
    MergeCount = amptklib.countfastq(merged_reads)
    amptklib.losslessTrim(merged_reads, FwdPrimer, RevPrimer, args.primer_mismatch, args.trim_len, args.pad, args.min_len, DemuxOut)
    FinalCount = amptklib.countfastq(DemuxOut)
    TooShort = MergeCount - FinalCount
    stats = os.path.join(tmpdir, base+'.stats')
    with open(stats, 'w') as counts:
    	counts.write("%i,%i,%i,%i,%i,%i\n" % (Total, BCFound, ForPrimerCount, RevPrimerCount, TooShort, FinalCount))

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: processReadsPE(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))
        
def main(args):
	global FwdPrimer, RevPrimer, Barcodes, tmpdir, usearch
	parser=argparse.ArgumentParser(prog='amptk-process_illumina_raw.py', 
		usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
		description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-f','--forward', dest='fastq', required=True, help='Illumina FASTQ R1 reads')
	parser.add_argument('-r', '--reverse', required=True, help='Illumina FASTQ R2 reads')
	parser.add_argument('-i', '--index', nargs='+', required=True, help='Illumina FASTQ index reads')
	parser.add_argument('-m', '--mapping_file', help='QIIME-like mapping file')
	parser.add_argument('--read_length', type=int, help='Read length, i.e. 2 x 300 bp = 300')
	parser.add_argument('-o','--out', dest="out", default='illumina_out', help='Base name for output')
	parser.add_argument('--fwd_primer', dest="F_primer", default='515FB', help='Forward Primer')
	parser.add_argument('--rev_primer', dest="R_primer", default='806RB', help='Reverse Primer')
	parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
	parser.add_argument('--barcode_mismatch', default=0, type=int, help='Number of mis-matches in barcode')
	parser.add_argument('--barcode_fasta', help='FASTA file containing Barcodes (Names & Sequences)')
	parser.add_argument('--rescue_forward', default='on', choices=['on', 'off'], help='Rescue Not-merged forward reads')
	parser.add_argument('--barcode_rev_comp', action='store_true', help='Reverse complement barcode sequences')
	parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
	parser.add_argument('-l','--trim_len', default=300, type=int, help='Trim length for reads')
	parser.add_argument('-p','--pad', default='off', choices=['on', 'off'], help='Pad with Ns to a set length')
	parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
	parser.add_argument('--cleanup', action='store_true', help='remove intermediate files')
	parser.add_argument('--merge_method', default='usearch', choices=['usearch', 'vsearch'], help='Software to use for PE read merging')
	args=parser.parse_args(args)

	args.out = re.sub(r'\W+', '', args.out)

	log_name = args.out+'.amptk-demux.log'
	if os.path.isfile(log_name):
		os.remove(log_name)

	amptklib.setupLogging(log_name)
	FNULL = open(os.devnull, 'w')
	cmd_args = " ".join(sys.argv)+'\n'
	amptklib.log.debug(cmd_args)
	print("-------------------------------------------------------")

	#initialize script, log system info and usearch version
	amptklib.SystemInfo()
	#get version of amptk
	usearch = args.usearch
	amptklib.versionDependencyChecks(usearch)

	#get number of CPUs to use
	if not args.cpus:
		cpus = multiprocessing.cpu_count()
	else:
		cpus = args.cpus

	#create tmpdir
	tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)

	#parse a mapping file or a barcode fasta file, primers, etc get setup
	#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
	barcode_file = args.out + ".barcodes_used.fa"
	if os.path.isfile(barcode_file):
		os.remove(barcode_file)

	#check if mapping file passed, use this if present, otherwise use command line arguments
	SampleData = {}
	Barcodes = {}
	RevBarcodes = {}
	FwdPrimer = ''
	RevPrimer = ''
	if args.mapping_file:
		if not os.path.isfile(args.mapping_file):
			amptklib.log.error("Mapping file not found: %s" % args.mapping_file)
			sys.exit(1)
		SampleData, Barcodes, RevBarcodes, FwdPrimer, RevPrimer = amptklib.parseMappingFileNEW(args.mapping_file)  
	else: #no mapping file, so create dictionaries from barcode fasta files
		if not args.barcode_fasta:
			amptklib.log.error("You did not specify a --barcode_fasta or --mapping_file, one is required")
			sys.exit(1)
		else:
			shutil.copyfile(args.barcode_fasta, barcode_file)
			Barcodes = amptklib.fasta2barcodes(barcode_file, False)                
		
	if FwdPrimer == '' or RevPrimer == '':
		#parse primers here so doesn't conflict with mapping primers
		#look up primer db otherwise default to entry
		if args.F_primer in amptklib.primer_db:
			FwdPrimer = amptklib.primer_db.get(args.F_primer)
			amptklib.log.info("{:} fwd primer found in AMPtk primer db, setting to: {:}".format(args.F_primer, FwdPrimer))
		else:
			FwdPrimer = args.F_primer
			amptklib.log.info("{:} fwd primer not found in AMPtk primer db, assuming it is actual primer sequence.".format(args.F_primer))
		if args.R_primer in amptklib.primer_db:
			RevPrimer = amptklib.primer_db.get(args.R_primer)
			amptklib.log.info("{:} rev primer found in AMPtk primer db, setting to: {:}".format(args.R_primer, RevPrimer))
		else:
			RevPrimer = args.R_primer
			amptklib.log.info("{:} rev primer not found in AMPtk primer db, assuming it is actual primer sequence.".format(args.R_primer))

	#if still no primers set, then exit
	if FwdPrimer == '' or RevPrimer == '':
		amptklib.log.error("Please provide primer sequences via --fwd_primer and --rev_primer")
		sys.exit(1)

	#if barcodes_rev_comp passed then reverse complement the keys in mapdict
	if args.barcode_rev_comp:
		amptklib.log.info("Reverse complementing barcode sequences")
		backupDict = Barcodes
		Barcodes = {}
		for k,v in list(backupDict.items()):
			RCkey = amptklib.RevComp(v)
			Barcodes[k] = RCkey

	amptklib.log.info("Loading %i samples from mapping file" % len(Barcodes))
	amptklib.log.info('FwdPrimer: {:}  RevPrimer: {:}'.format(FwdPrimer, RevPrimer))
	amptklib.log.info('Dropping reads less than {:} bp and setting lossless trimming to {:} bp.'.format(args.min_len, args.trim_len))

	#rename reads according to indexes
	if not amptklib.PEandIndexCheck(args.fastq, args.reverse, args.index[0]): #check they are all same length
		amptklib.log.error("FASTQ input malformed, read numbers do not match")
		sys.exit(1)
	amptklib.log.info("Loading FASTQ Records")
	NumSeqs = amptklib.countfastq(args.fastq)
	if cpus > 1:
		amptklib.log.info("Splitting FASTQ files over {:} cpus".format(cpus))
		amptklib.split_fastqPEandI(args.fastq, args.reverse, args.index[0], NumSeqs, tmpdir, cpus*2)
		file_list = []
		for file in os.listdir(tmpdir):
			if file.endswith('.fq'):
				filepart = os.path.join(tmpdir, file.split('_R')[0])
				if not filepart in file_list:
					file_list.append(filepart)

		amptklib.log.info("Mapping indexes to reads and renaming PE reads")
		amptklib.runMultiProgress(safe_run, file_list, cpus, args=args)
	else:
		amptklib.log.info("Mapping indexes to reads and renaming PE reads")
		shutil.copyfile(args.fastq, os.path.join(tmpdir, 'chunk_R1.fq'))
		shutil.copyfile(args.reverse, os.path.join(tmpdir, 'chunk_R2.fq'))
		shutil.copyfile(args.index[0], os.path.join(tmpdir, 'chunk_R3.fq'))
		processReadsPE(os.path.join(tmpdir, 'chunk'), args=args)

	print("-------------------------------------------------------")
	#Now concatenate all of the demuxed files together
	amptklib.log.info("Concatenating Demuxed Files")

	tmpDemux = os.path.join(tmpdir, args.out + '.demux.fq')
	with open(tmpDemux, 'wb') as outfile:
		for filename in glob.glob(os.path.join(tmpdir,'*.demux.fq')):
			if filename == tmpDemux:
				continue
			with open(filename, 'r') as readfile:
				shutil.copyfileobj(readfile, outfile)
	#parse the stats
	finalstats = [0,0,0,0,0,0]
	for file in os.listdir(tmpdir):
		if file.endswith('.stats'):
			with open(os.path.join(tmpdir, file), 'r') as statsfile:
				line = statsfile.readline()
				line = line.replace('\n', '')
				newstats = line.split(',')
				newstats = [int(i) for i in newstats]
				for x, num in enumerate(newstats):
					finalstats[x] += num

	#finally reindex output
	#last thing is to re-number of reads as it is possible they could have same name from multitprocessor split
	Demux = args.out + '.demux.fq'
	amptklib.fastqreindex(tmpDemux, Demux)
	amptklib.SafeRemove(tmpDemux)

	#output stats of the run
	amptklib.log.info('{0:,}'.format(finalstats[0])+' total reads')
	amptklib.log.info('{0:,}'.format(finalstats[0] - finalstats[1])+' discarded no index match')
	amptklib.log.info('{0:,}'.format(finalstats[2])+' Fwd Primer found, {0:,}'.format(finalstats[3])+ ' Rev Primer found')
	amptklib.log.info('{0:,}'.format(finalstats[4])+' discarded too short (< %i bp)' % args.min_len)
	amptklib.log.info('{0:,}'.format(finalstats[5])+' valid output reads')

	#now loop through data and find barcoded samples, counting each.....
	BarcodeCount = {}
	with open(Demux, 'r') as input:
		header = itertools.islice(input, 0, None, 4)
		for line in header:
			ID = line.split("=",1)[-1].split(";")[0]
			if ID not in BarcodeCount:
				BarcodeCount[ID] = 1
			else:
				BarcodeCount[ID] += 1

	#now let's count the barcodes found and count the number of times they are found.
	barcode_counts = "%30s:  %s" % ('Sample', 'Count')
	for k,v in natsorted(list(BarcodeCount.items()), key=lambda k_v: k_v[1], reverse=True):
		barcode_counts += "\n%30s:  %s" % (k, str(BarcodeCount[k]))
	amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

	#create mapping file if one doesn't exist
	genericmapfile = args.out + '.mapping_file.txt'
	amptklib.CreateGenericMappingFile(Barcodes, {}, FwdPrimer, RevPrimer, genericmapfile, BarcodeCount)

	#compress the output to save space
	FinalDemux = Demux+'.gz'
	amptklib.Fzip(Demux, FinalDemux, cpus)
	amptklib.removefile(Demux)

	if args.cleanup:
		amptklib.SafeRemove(tmpdir)

	#get file size
	filesize = os.path.getsize(FinalDemux)
	readablesize = amptklib.convertSize(filesize)
	amptklib.log.info("Output file:  %s (%s)" % (FinalDemux, readablesize))
	amptklib.log.info("Mapping file: %s" % genericmapfile)
	print("-------------------------------------------------------")
	if 'darwin' in sys.platform:
		print(col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (FinalDemux))
	else:
		print("\nExample of next cmd: amptk cluster -i %s -o out\n" % (FinalDemux))   

if __name__ == "__main__":
	main(args)