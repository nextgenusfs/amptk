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
import time
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

def processRead(input, args=False):
	base = os.path.basename(input).split('.')[0]
	''''
	FwdPrimer = dictOptions['FwdPrimer']
	RevPrimer = dictOptions['RevPrimer']
	Barcodes = dictOptions['Barcodes']
	args = dictOptions['args']
	'''
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
				revalign = edlib.align(RevPrimer, Seq, mode="HW", task="locations", k=args.primer_mismatch, additionalEqualities=amptklib.degenNuc)
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
	global FwdPrimer, RevPrimer, Barcodes, tmpdir
	parser=argparse.ArgumentParser(prog='amptk-process_ion.py', usage="%(prog)s [options] -i file.fastq\n%(prog)s -h for help menu",
		description='''Script finds barcodes, strips forward and reverse primers, relabels, and then trim/pads reads to a set length''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--fastq','--sff', '--fasta', '--bam', dest='fastq', required=True, help='BAM/FASTQ/SFF/FASTA file')
	parser.add_argument('-q','--qual', help='QUAL file (if -i is FASTA)')
	parser.add_argument('-o','--out', dest="out", default='ion', help='Base name for output')
	parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7-ion', help='Forward Primer')
	parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer')
	parser.add_argument('-m','--mapping_file', help='Mapping file: QIIME format can have extra meta data columns')
	parser.add_argument('-p','--pad', default='off', choices=['on', 'off'], help='Pad with Ns to a set length')
	parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
	parser.add_argument('--barcode_mismatch', default=0, type=int, help='Number of mis-matches in barcode')
	parser.add_argument('--barcode_fasta', default='ionxpress', help='FASTA file containing Barcodes (Names & Sequences)')
	parser.add_argument('--reverse_barcode', help='FASTA file containing 3 prime Barocdes')
	parser.add_argument('-b','--list_barcodes', dest="barcodes", default='all', help='Enter Barcodes used separated by commas')
	parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
	parser.add_argument('-l','--trim_len', default=300, type=int, help='Trim length for reads')
	parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
	parser.add_argument('--mult_samples', dest="multi", default='False', help='Combine multiple samples (i.e. FACE1)')
	parser.add_argument('--ion', action='store_true', help='Input data is Ion Torrent')
	parser.add_argument('--454', action='store_true', help='Input data is 454')
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
	if args.mapping_file:
		if not os.path.isfile(args.mapping_file):
			amptklib.log.error("Mapping file not found: %s" % args.mapping_file)
			sys.exit(1)
		SampleData, Barcodes, RevBarcodes, FwdPrimer, RevPrimer = amptklib.parseMappingFileNEW(args.mapping_file)  
		genericmapfile = args.mapping_file
	else: #no mapping file, so create dictionaries from barcode fasta files
		if args.barcode_fasta == 'ionxpress':
			#get script path and barcode file name
			pgm_barcodes = os.path.join(os.path.dirname(amptklib.__file__), 'DB', 'ionxpress_barcodes.fa')
		elif args.barcode_fasta == 'ioncode':
			pgm_barcodes = os.path.join(os.path.dirname(amptklib.__file__), 'DB', 'ioncode_barcodes.fa')
		if args.barcode_fasta == 'ionxpress' or args.barcode_fasta == 'ioncode':
			if args.barcodes == "all":
				if args.multi == 'False':
					shutil.copyfile(pgm_barcodes, barcode_file)
				else:
					with open(barcode_file, 'w') as barcodeout:
						with open(pgm_barcodes, 'r') as input:
							for rec in SeqIO.parse(input, 'fasta'):
								outname = args.multi+'.'+rec.id
								barcodeout.write(">%s\n%s\n" % (outname, rec.seq))
			else:
				bc_list = args.barcodes.split(",")
				inputSeqFile = open(pgm_barcodes, "rU")
				SeqRecords = SeqIO.to_dict(SeqIO.parse(inputSeqFile, "fasta"))
				for rec in bc_list:
					name = "BC." + rec
					seq = SeqRecords[name].seq
					if args.multi != 'False':
						outname = args.multi+'.'+name
					else:
						outname = name
					outputSeqFile = open(barcode_file, "a")
					outputSeqFile.write(">%s\n%s\n" % (outname, seq))
				outputSeqFile.close()
				inputSeqFile.close()
		else:
			#check for multi_samples and add if necessary
			if args.multi == 'False':
				shutil.copyfile(args.barcode_fasta, barcode_file)
				if args.reverse_barcode:
					shutil.copyfile(args.reverse_barcode,rev_barcode_file)
			else:
				with open(barcode_file, 'w') as barcodeout:
					with open(args.barcode_fasta, 'r') as input:
						for rec in SeqIO.parse(input, 'fasta'):
							outname = args.multi+'.'+rec.id
							barcodeout.write(">%s\n%s\n" % (outname, rec.seq))
				if args.reverse_barcode:
					with open(rev_barcode_file, 'w') as barcodeout:
						with open(args.reverse_barcode, 'r') as input:
							for rec in SeqIO.parse(input, 'fasta'):
								outname = args.multi+'.'+rec.id
								barcodeout.write(">%s\n%s\n" % (outname, rec.seq))                   
	
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
	if gzip_list:
		amptklib.log.info("Gzipped input files detected, uncompressing")
		for file in gzip_list:
			file_out = file.replace('.gz', '')
			amptklib.Funzip(file, file_out, cpus)
		args.fastq = args.fastq.replace('.gz', '')
	 
	#if SFF file passed, convert to FASTQ with biopython
	if args.fastq.endswith(".sff"):
		if args.barcode_fasta == 'ionxpress':
			if not args.mapping_file:
				amptklib.log.error("You did not specify a --barcode_fasta or --mapping_file, one is required for 454 data")
				sys.exit(1)
		amptklib.log.info("SFF input detected, converting to FASTQ")
		SeqIn = args.out + '.sff.extract.fastq'
		SeqIO.convert(args.fastq, "sff-trim", SeqIn, "fastq")
	elif args.fastq.endswith(".fas") or args.fastq.endswith(".fasta") or args.fastq.endswith(".fa"):
		if not args.qual:
			amptklib.log.error("FASTA input detected, however no QUAL file was given.  You must have FASTA + QUAL files")
			sys.exit(1)
		else:
			if args.barcode_fasta == 'ionxpress':
				if not args.mapping_file:
					amptklib.log.error("You did not specify a --barcode_fasta or --mapping_file, one is required for 454 data")
					sys.exit(1)
			SeqIn = args.out + '.fastq'
			amptklib.log.info("FASTA + QUAL detected, converting to FASTQ")
			amptklib.faqual2fastq(args.fastq, args.qual, SeqIn)
	elif args.fastq.endswith('.bam'):
		#so we can convert natively with pybam, however it is 10X slower than bedtools/samtools
		#since samtools is fastest, lets use that if exists, if not then bedtools, else default to pybam
		amptklib.log.info("Converting Ion Torrent BAM file to FASTQ")
		SeqIn = args.out+'.fastq'
		if amptklib.which('samtools'):
			cmd = ['samtools', 'fastq', '-@', str(cpus), args.fastq]
			amptklib.runSubprocess2(cmd, amptklib.log, SeqIn)
		else:
			if amptklib.which('bedtools'):
				cmd = ['bedtools', 'bamtofastq', '-i', args.fastq, '-fq', SeqIn]
				amptklib.runSubprocess(cmd, amptklib.log)
			else: #default to pybam
				amptklib.bam2fastq(args.fastq, SeqIn)
	else:        
		SeqIn = args.fastq

	#start here to process the reads, first reverse complement the reverse primer
	catDemux = args.out + '.demux.fq'
	origRevPrimer = RevPrimer
	RevPrimer = amptklib.RevComp(RevPrimer)
	amptklib.log.info("Foward primer: %s,  Rev comp'd rev primer: %s" % (FwdPrimer, RevPrimer))

	#then setup barcode dictionary
	if len(Barcodes) < 1:
		Barcodes = amptklib.fasta2barcodes(barcode_file, False)

	#setup for looking for reverse barcode
	if len(RevBarcodes) < 1 and args.reverse_barcode:
		if not os.path.isfile(args.reverse_barcode):
			amptklib.log.info("Reverse barcode is not a valid file, exiting")
			sys.exit(1) 
		shutil.copyfile(args.reverse_barcode, rev_barcode_file)
		RevBarcodes = amptklib.fasta2barcodes(rev_barcode_file, True)

	#Count FASTQ records
	amptklib.log.info("Loading FASTQ Records")
	orig_total = amptklib.countfastq(SeqIn)
	size = amptklib.checkfastqsize(SeqIn)
	readablesize = amptklib.convertSize(size)
	amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

	#create tmpdir and split input into n cpus
	tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)

	amptklib.log.info('Dropping reads less than {:} bp and setting lossless trimming to {:} bp.'.format(args.min_len, args.trim_len))

	if cpus > 1:
		#split fastq file
		amptklib.log.info("Splitting FASTQ files over {:} cpus".format(cpus))
		amptklib.split_fastq(SeqIn, orig_total, tmpdir, cpus*2)    
		#now get file list from tmp folder
		file_list = []
		for file in os.listdir(tmpdir):
			if file.endswith(".fq"):
				file = os.path.join(tmpdir, file)
				file_list.append(file)
		#finally process reads over number of cpus
		amptklib.runMultiProgress(processRead, file_list, cpus, args=args)
	else:
		shutil.copyfile(SeqIn, os.path.join(tmpdir, 'chunk.fq'))
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
		
	#clean up tmp folder
	shutil.rmtree(tmpdir)

	#last thing is to re-number of reads as it is possible they could have same name from multitprocessor split
	amptklib.fastqreindex(tmpDemux, catDemux)
	os.remove(tmpDemux)
	
	amptklib.log.info('{0:,}'.format(finalstats[0])+' total reads')
	if args.reverse_barcode:
		amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1]-finalstats[2]-finalstats[4])+' valid Fwd and Rev Barcodes')
	else:
		amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1])+' valid Barcode')
		amptklib.log.info('{0:,}'.format(finalstats[0]-finalstats[1]-finalstats[2])+' Fwd Primer found, {0:,}'.format(finalstats[3])+ ' Rev Primer found')
	amptklib.log.info('{0:,}'.format(finalstats[5])+' discarded too short (< %i bp)' % args.min_len)
	amptklib.log.info('{0:,}'.format(finalstats[6])+' valid output reads')

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

	#create a generic mappingfile for downstream processes
	genericmapfile = args.out + '.mapping_file.txt'
	if not args.mapping_file:
		amptklib.CreateGenericMappingFile(Barcodes, RevBarcodes, FwdPrimer, origRevPrimer, genericmapfile, BarcodeCount)
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