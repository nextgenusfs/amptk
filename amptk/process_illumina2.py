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
    orientR1 = os.path.join(tmpdir, base+'_R1.oriented.fq')
    orientR2 = os.path.join(tmpdir, base+'_R2.oriented.fq')
    trim_forward = os.path.join(tmpdir, base+'_R1.trimmed.fq')
    trim_reverse = os.path.join(tmpdir, base+'_R2.trimmed.fq')
    merged_reads = os.path.join(tmpdir, base+'.merged.fq')
    DemuxOut = os.path.join(tmpdir, base+'.demux.fq')
    StatsOut = os.path.join(tmpdir, base+'.stats')
    RL = amptklib.GuessRL(forward_reads)
    Total, Correct, Flip, Drop = amptklib.illuminaReorient(forward_reads, reverse_reads, FwdPrimer, RevPrimer, args.primer_mismatch, RL, orientR1, orientR2)
    amptklib.log.debug('Re-oriented PE reads for {:}: {:,} total, {:,} correct, {:,} flipped, {:,} dropped.'.format(base, Total, Correct, Flip, Drop))
    if args.barcode_not_anchored:
        amptklib.demuxIlluminaPE2(orientR1, orientR2, FwdPrimer, RevPrimer, SampleData, Barcodes, RevBarcodes, args.barcode_mismatch, args.primer_mismatch, trim_forward, trim_reverse, StatsOut)
    else:
        amptklib.demuxIlluminaPE(orientR1, orientR2, FwdPrimer, RevPrimer, SampleData, Barcodes, RevBarcodes, args.barcode_mismatch, args.primer_mismatch, trim_forward, trim_reverse, StatsOut)
    if args.full_length:
        amptklib.MergeReadsSimple(trim_forward, trim_reverse, '.', DemuxOut, args.min_len, usearch, 'off', args.merge_method)
    else:
        amptklib.MergeReadsSimple(trim_forward, trim_reverse, '.', merged_reads, args.min_len, usearch, 'on', args.merge_method)
        amptklib.losslessTrim(merged_reads, FwdPrimer, RevPrimer, args.primer_mismatch, args.trim_len, args.pad, args.min_len, DemuxOut) 
    amptklib.SafeRemove(orientR1)
    amptklib.SafeRemove(orientR2)
    amptklib.SafeRemove(forward_reads)
    amptklib.SafeRemove(reverse_reads)
    amptklib.SafeRemove(merged_reads)

def processRead(input, args=False):
    base = os.path.basename(input).split('.')[0]
    PL = len(FwdPrimer)
    RL = len(RevPrimer)
    DemuxOut = os.path.join(tmpdir, base+'.demux.fq')
    StatsOut = os.path.join(tmpdir, base+'.stats')
    Total = 0
    NoBarcode = 0
    NoRevBarcode = 0
    NoPrimer = 0
    TooShort = 0
    RevPrimerFound = 0
    ValidSeqs = 0
    RCrevPrimer = amptklib.RevComp(RevPrimer)
    with open(StatsOut, 'w') as counts:
        with open(DemuxOut, 'w') as out:   
            for title, seq, qual in FastqGeneralIterator(open(input)):
                Total += 1
                #look for barcode, trim it off
                Barcode, BarcodeLabel = amptklib.AlignBarcode(seq, Barcodes, args.barcode_mismatch)
                if Barcode == "":
                    NoBarcode += 1
                    continue
                BarcodeLength = len(Barcode)
                Seq = seq[BarcodeLength:]
                Qual = qual[BarcodeLength:]
                #now search for forward primer
                foralign = edlib.align(FwdPrimer, Seq, mode="HW", k=args.primer_mismatch, additionalEqualities=amptklib.degenNuc)
                if foralign["editDistance"] < 0:
                    NoPrimer += 1
                    continue
                ForTrim = foralign["locations"][0][1]+1   
                #now search for reverse primer
                revalign = edlib.align(RCrevPrimer, Seq, mode="HW", task="locations", k=args.primer_mismatch, additionalEqualities=amptklib.degenNuc)
                if revalign["editDistance"] >= 0:  #reverse primer was found
                    RevPrimerFound += 1 
                    #location to trim sequences
                    RevTrim = revalign["locations"][0][0]                
                    #determine reverse barcode
                    if args.reverse_barcode:
                        RevBCdiffs = 0
                        BCcut = revalign["locations"][0][1]
                        CutSeq = Seq[BCcut:]
                        RevBarcode, RevBarcodeLabel = amptklib.AlignRevBarcode(CutSeq, RevBarcodes, args.barcode_mismatch)
                        if RevBarcode == "":
                            NoRevBarcode += 1
                            continue
                        BarcodeLabel = BarcodeLabel+':-:'+RevBarcodeLabel                       
                    #now trim record remove forward and reverse reads
                    Seq = Seq[ForTrim:RevTrim]
                    Qual = Qual[ForTrim:RevTrim]
                    #since found reverse primer, now also need to pad/trim
                    if not args.full_length:
                        #check minimum length here or primer dimer type sequences will get padded with Ns
                        if len(Seq) < int(args.min_len):
                            TooShort += 1
                            continue
                        if len(Seq) < args.trim_len and args.pad == 'on':
                            pad = args.trim_len - len(Seq)
                            Seq = Seq + pad*'N'
                            Qual = Qual +pad*'I'
                        else: #len(Seq) > args.trim_len:
                            Seq = Seq[:args.trim_len]
                            Qual = Qual[:args.trim_len]
                else:
                    #trim record, did not find reverse primer
                    if args.full_length: #if full length then move to next record
                        continue
                    #trim away forward primer
                    Seq = Seq[ForTrim:]
                    Qual = Qual[ForTrim:]
                    #check length and trim, throw away if too short as it was bad read
                    if len(Seq) < args.trim_len:
                        TooShort += 1
                        continue
                    Seq = Seq[:args.trim_len]
                    Qual = Qual[:args.trim_len]
                #check minimum length
                if len(Seq) < int(args.min_len):
                    TooShort += 1
                    continue
                ValidSeqs += 1
                #rename header
                Name = 'R_'+str(ValidSeqs)+';barcodelabel='+BarcodeLabel+';'
                out.write("@%s\n%s\n+\n%s\n" % (Name, Seq, Qual))
            counts.write('%i,%i,%i,%i,%i,%i,%i\n' % (Total, NoBarcode, NoPrimer, RevPrimerFound, NoRevBarcode, TooShort, ValidSeqs))    
    
def main(args):
	global FwdPrimer, RevPrimer, SampleData, Barcodes, RevBarcodes, tmpdir, usearch
	parser=argparse.ArgumentParser(prog='amptk-process_ion.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
		description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--fastq', dest='fastq', required=True, help='FASTQ R1 file')
	parser.add_argument('--reverse', help='Illumina R2 reverse reads')
	parser.add_argument('-o','--out', dest="out", default='illumina2', help='Base name for output')
	parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer')
	parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer')
	parser.add_argument('-m','--mapping_file', help='Mapping file: QIIME format can have extra meta data columns')
	parser.add_argument('-p','--pad', default='off', choices=['on', 'off'], help='Pad with Ns to a set length')
	parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
	parser.add_argument('--barcode_mismatch', default=0, type=int, help='Number of mis-matches in barcode')
	parser.add_argument('--barcode_fasta', help='FASTA file containing Barcodes (Names & Sequences)')
	parser.add_argument('--barcode_not_anchored', action='store_true', help='Barcodes (indexes) are not at start of reads')
	parser.add_argument('--reverse_barcode', help='FASTA file containing 3 prime Barocdes')
	parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
	parser.add_argument('-l','--trim_len', default=300, type=int, help='Trim length for reads')
	parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
	parser.add_argument('--merge_method', default='usearch', choices=['usearch', 'vsearch'], help='Software to use for PE read merging')
	parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH EXE')
	args=parser.parse_args(args)    

	args.out = re.sub(r'\W+', '', args.out)

	log_name = args.out + '.amptk-demux.log'
	if os.path.isfile(log_name):
		os.remove(log_name)
	FNULL = open(os.devnull, 'w')
	amptklib.setupLogging(log_name)
	cmd_args = " ".join(sys.argv)+'\n'
	amptklib.log.debug(cmd_args)
	print("-------------------------------------------------------")

	#initialize script, log system info and usearch version
	amptklib.SystemInfo()
	#Do a version check
	usearch = args.usearch
	amptklib.versionDependencyChecks(usearch)

	#get number of CPUs to use
	if not args.cpus:
		cpus = multiprocessing.cpu_count()
	else:
		cpus = args.cpus

	#parse a mapping file or a barcode fasta file, primers, etc get setup
	#dealing with Barcodes, get ion barcodes or parse the barcode_fasta argument
	barcode_file = args.out + ".barcodes_used.fa"
	rev_barcode_file = args.out + '.revbarcodes_used.fa'
	amptklib.SafeRemove(barcode_file)
	amptklib.SafeRemove(rev_barcode_file)

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
			if args.reverse_barcode:
				shutil.copyfile(args.reverse_barcode, rev_barcode_file)
				RevBarcodes = amptklib.fasta2barcodes(rev_barcode_file, False)                   
	
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

	#check if input is compressed
	gzip_list = []
	if args.fastq.endswith('.gz'):
		gzip_list.append(os.path.abspath(args.fastq))
	if args.reverse:
		if args.reverse.endswith('.gz'):
			gzip_list.append(os.path.abspath(args.reverse))
	if gzip_list:
		amptklib.log.info("Gzipped input files detected, uncompressing")
		for file in gzip_list:
			file_out = file.replace('.gz', '')
			amptklib.Funzip(file, file_out, cpus)
		args.fastq = args.fastq.replace('.gz', '')
		if args.reverse:
			args.reverse = args.reverse.replace('.gz', '')

	#Count FASTQ records
	amptklib.log.info("Loading FASTQ Records")
	orig_total = amptklib.countfastq(args.fastq)
	size = amptklib.checkfastqsize(args.fastq)
	readablesize = amptklib.convertSize(size*2)
	amptklib.log.info('{:,} reads ({:})'.format(orig_total, readablesize))

	#output barcodes/samples
	amptklib.log.info('Searching for {:} forward barcodes and {:} reverse barcodes'.format(len(Barcodes), len(RevBarcodes)))

	#create tmpdir and split input into n cpus
	tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)
	
	#tell user about number of cores using
	amptklib.log.info("Splitting FASTQ files over {:} cpus".format(cpus))

	if args.reverse:
		amptklib.log.info("Demuxing PE Illumina reads; FwdPrimer: {:} RevPrimer: {:}".format(FwdPrimer, RevPrimer))
	else:
		amptklib.log.info("Demuxing SE Illumina reads; FwdPrimer: {:} RevPrimer: {:}".format(FwdPrimer, amptklib.RevComp(RevPrimer)))

	amptklib.log.info('Dropping reads less than {:} bp and setting lossless trimming to {:} bp.'.format(args.min_len, args.trim_len))

	if cpus > 1:
		if args.reverse:
			amptklib.split_fastqPE(args.fastq, args.reverse, orig_total, tmpdir, cpus*4)
			file_list = []
			for file in os.listdir(tmpdir):
				if file.endswith('.fq'):
					filepart = os.path.join(tmpdir, file.split('_R')[0])
					if not filepart in file_list:
						file_list.append(filepart)
			amptklib.runMultiProgress(processReadsPE, file_list, cpus, args=args)               
		else:
			#split fastq file
			amptklib.split_fastq(args.fastq, orig_total, tmpdir, cpus*4)    
			#now get file list from tmp folder
			file_list = []
			for file in os.listdir(tmpdir):
				if file.endswith(".fq"):
					file = os.path.join(tmpdir, file)
					file_list.append(file)
			#finally process reads over number of cpus
			amptklib.runMultiProgress(processRead, file_list, cpus, args=args)
	else:
		if args.reverse:
			shutil.copyfile(args.fastq, os.path.join(tmpdir, 'chunk_R1.fq'))
			shutil.copyfile(args.reverse, os.path.join(tmpdir, 'chunk_R2.fq'))
			processReadsPE(os.path.join(tmpdir, 'chunk'), args=args)
		else:
			shutil.copyfile(args.fastq, os.path.join(tmpdir, 'chunk.fq'))
			processRead(os.path.join(tmpdir, 'chunk.fq'), args=args)

	print("-------------------------------------------------------")
	#Now concatenate all of the demuxed files together
	amptklib.log.info("Concatenating Demuxed Files")

	tmpDemux = args.out + '.tmp.demux.fq'
	with open(tmpDemux, 'w') as outfile:
		for filename in glob.glob(os.path.join(tmpdir,'*.demux.fq')):
			if filename == tmpDemux:
				continue
			with open(filename, 'r') as readfile:
				shutil.copyfileobj(readfile, outfile)
	if args.reverse:
		#parse the stats
		finalstats = [0,0,0,0,0,0]
		for file in os.listdir(tmpdir):
			if file.endswith('.stats'):
				with open(os.path.join(tmpdir, file), 'r') as statsfile:
					line = statsfile.readline()
					line = line.rstrip()
					newstats = line.split(',')
					newstats = [int(i) for i in newstats]
					for x, num in enumerate(newstats):
						finalstats[x] += num
	
		amptklib.log.info('{0:,}'.format(finalstats[0])+' total reads')
		amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1]-finalstats[3])+' valid Barcodes')
		amptklib.log.info('{0:,}'.format(finalstats[5])+' valid output reads (Barcodes and Primers)')
	else:
		#parse the stats
		finalstats = [0,0,0,0,0,0,0]
		for file in os.listdir(tmpdir):
			if file.endswith('.stats'):
				with open(os.path.join(tmpdir, file), 'r') as statsfile:
					line = statsfile.readline()
					line = line.rstrip()
					newstats = line.split(',')
					newstats = [int(i) for i in newstats]
					for x, num in enumerate(newstats):
						finalstats[x] += num
			
		amptklib.log.info('{0:,}'.format(finalstats[0])+' total reads')
		if args.reverse_barcode:
			amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1]-finalstats[2]-finalstats[4])+' valid Fwd and Rev Barcodes')
		else:
			amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1])+' valid Barcode')
			amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1]-finalstats[2])+' Fwd Primer found, {0:,}'.format(finalstats[3])+ ' Rev Primer found')
		amptklib.log.info('{0:,}'.format(finalstats[5])+' discarded too short (< %i bp)' % args.min_len)
		amptklib.log.info('{0:,}'.format(finalstats[6])+' valid output reads')


	#clean up tmp folder
	amptklib.SafeRemove(tmpdir)

	#last thing is to re-number of reads as it is possible they could have same name from multitprocessor split
	catDemux = args.out+'.demux.fq'
	amptklib.fastqreindex(tmpDemux, catDemux)
	amptklib.SafeRemove(tmpDemux)
	#now loop through data and find barcoded samples, counting each.....
	BarcodeCount = {}
	with open(catDemux, 'r') as input:
		header = itertools.islice(input, 0, None, 4)
		for line in header:
			ID = line.split("=",1)[-1].split(";")[0]
			if ID not in BarcodeCount:
				BarcodeCount[ID] = 1
			else:
				BarcodeCount[ID] += 1

	#now let's count the barcodes found and count the number of times they are found.
	barcode_counts = "%22s:  %s" % ('Sample', 'Count')
	for k,v in natsorted(list(BarcodeCount.items()), key=lambda k_v: k_v[1], reverse=True):
		barcode_counts += "\n%22s:  %s" % (k, str(BarcodeCount[k]))
	amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

	genericmapfile = args.out + '.mapping_file.txt'
	if not args.mapping_file:
		#create a generic mappingfile for downstream processes
		amptklib.CreateGenericMappingFile(Barcodes, RevBarcodes, FwdPrimer, RevPrimer, genericmapfile, BarcodeCount)
	else:
		amptklib.updateMappingFile(args.mapping_file, BarcodeCount, genericmapfile)
	#compress the output to save space
	FinalDemux = catDemux+'.gz'
	amptklib.Fzip(catDemux, FinalDemux, cpus)
	amptklib.removefile(catDemux)
	if gzip_list:
		for file in gzip_list:
			file = file.replace('.gz', '')
			amptklib.removefile(file)

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
