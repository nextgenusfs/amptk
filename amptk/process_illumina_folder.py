#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os
import sys
import argparse
import shutil
import subprocess
import glob
import math
import logging
import gzip
import inspect
import multiprocessing
import itertools
import re
from natsort import natsorted
import edlib
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from amptk import amptklib


def processSEreads(input, args=False):
    #input is expected to be a FASTQ file
    #local variables that need to be previously declared: FwdPrimer, RevPrimer
    inputPath = os.path.abspath(input)
    Name = os.path.basename(inputPath).split(".fq",-1)[0]
    DemuxOut = os.path.join(args.out, Name + '.demux.fq')
    Sample = Name.split('_')[0]
    StatsOut = os.path.join(args.out, Name+'.stats')
    Total = 0
    NoPrimer = 0
    TooShort = 0
    RevPrimerFound = 0
    ValidSeqs = 0
    multiHits = 0
    RCRevPrimer = amptklib.RevComp(RevPrimer)
    with open(StatsOut, 'w') as counts:
        with open(DemuxOut, 'w') as out:
            for title, seq, qual in FastqGeneralIterator(open(inputPath)):
                Total += 1
                #first thing is look for forward primer, if found trim it off
                foralign = edlib.align(FwdPrimer, seq, mode="HW", k=args.primer_mismatch, additionalEqualities=amptklib.degenNuc)
                #if require primer is on make finding primer in amplicon required if amplicon is larger than read length
                #if less than read length, can't enforce primer because could have been trimmed via staggered trim in fastq_mergepairs
                if len(foralign['locations']) > 1:
                    multiHits += 1
                    continue
                if args.primer == 'on':
                    if foralign["editDistance"] < 0:
                        NoPrimer += 1
                        continue
                    ForCutPos = foralign["locations"][0][1]+1
                    Seq = seq[ForCutPos:]
                    Qual = qual[ForCutPos:]
                else:
                    if foralign["editDistance"] >= 0:
                        ForCutPos = foralign["locations"][0][1]+1
                        Seq = seq[ForCutPos:]
                        Qual = qual[ForCutPos:]
                    else:
                        NoPrimer += 1
                        Seq = seq
                        Qual = qual
                #now look for reverse primer
                revalign = edlib.align(RCRevPrimer, Seq, mode="HW", task="locations", k=args.primer_mismatch, additionalEqualities=amptklib.degenNuc)
                if revalign["editDistance"] >= 0:
                    RevPrimerFound += 1
                    RevCutPos = revalign["locations"][0][0]
                    #location to trim sequences, trim seqs
                    Seq = Seq[:RevCutPos]
                    Qual = Qual[:RevCutPos]
                else:
                    if args.full_length:
                        continue
                #if full_length is passed, then only trim primers
                if not args.full_length:
                    #got here if primers were found they were trimmed
                    #now check seq length, pad if too short, trim if too long
                    if len(Seq) < args.min_len: #need this check here or primer dimers will get through
                        TooShort += 1
                        continue
                    if len(Seq) < args.trim_len and args.pad == 'on':
                        pad = args.trim_len - len(Seq)
                        Seq = Seq + pad*'N'
                        Qual = Qual + pad*'I'
                    else: #len(Seq) > args.trim_len:
                        Seq = Seq[:args.trim_len]
                        Qual = Qual[:args.trim_len]
                #got here, reads are primers trimmed and trim/padded, check length
                if len(Seq) < args.min_len:
                    TooShort += 1
                    continue
                ValidSeqs += 1     
                #now fix header
                Title = 'R_'+str(ValidSeqs)+';barcodelabel='+Sample+';'
                #now write to file
                out.write("@%s\n%s\n+\n%s\n" % (Title, Seq, Qual))
            ForPrimerFound = Total - NoPrimer - multiHits
            counts.write("%i,%i,%i,%i,%i,%i\n" % (Total, ForPrimerFound, RevPrimerFound, multiHits, TooShort, ValidSeqs))

def processPEreads(input, args=False):
    '''
    function for multiprocessing of the data, so take file list as input, need global forward/reverse list available
    '''
    for_reads, rev_reads = input
    if '_' in os.path.basename(for_reads):
        name = os.path.basename(for_reads).split("_")[0]
    else:
        name = os.path.basename(for_reads)
    amptklib.log.debug('{:}: {:} {:}'.format(name, for_reads, rev_reads))
    #for_reads = os.path.join(args.input, forwardRead)
    #rev_reads = os.path.join(args.input, reverseRead)
    StatsOut = os.path.join(args.out, name+'.stats')
    #if read length explicity passed use it otherwise measure it
    if args.read_length:
        read_length = args.read_length
    else:
        read_length = amptklib.GuessRL(for_reads)
    trimR1 = os.path.join(args.out, name+'_R1.fq')
    trimR2 = os.path.join(args.out, name+'_R2.fq')
    mergedReads = os.path.join(args.out, name+'.merged.fq')
    demuxReads = os.path.join(args.out, name+'.demux.fq')
    TotalCount, Written, DropMulti, FindForPrimer, FindRevPrimer = amptklib.stripPrimersPE(for_reads, rev_reads, read_length, name, FwdPrimer, RevPrimer, args.primer_mismatch, args.primer, args.full_length, trimR1, trimR2)
    MergedCount, PhixCleanedCount = amptklib.MergeReadsSimple(trimR1, trimR2, args.out, name+'.merged.fq', args.min_len, usearch, args.rescue_forward, args.merge_method)
    amptklib.losslessTrim(mergedReads, FwdPrimer, RevPrimer, args.primer_mismatch, args.trim_len, args.pad, args.min_len, demuxReads)
    FinalCount = amptklib.countfastq(demuxReads)
    TooShort = PhixCleanedCount - FinalCount
    with open(StatsOut, 'w') as counts:
        counts.write("%i,%i,%i,%i,%i,%i\n" % (TotalCount, FindForPrimer, FindRevPrimer, DropMulti, TooShort, FinalCount))

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: processPEreads(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def safe_run2(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: processSEreads(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))
        
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
	def __init__(self,prog):
		super(MyFormatter,self).__init__(prog,max_help_position=48)
class col(object):
	GRN = '\033[92m'
	END = '\033[0m'
	WARN = '\033[93m'

def main(args):
	global FwdPrimer, RevPrimer, usearch
	parser=argparse.ArgumentParser(prog='amptk-process_illumina_folder.py', usage="%(prog)s [options] -i folder",
		description='''Script that takes De-mulitplexed Illumina data from a folder and processes it for amptk (merge PE reads, strip primers, trim/pad to set length.''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--input', dest='input', required=True, help='Folder of Illumina Data')
	parser.add_argument('-o','--out', dest="out", default='amptk-illumina', help='Name for output folder')
	parser.add_argument('-m','--mapping_file', help='Mapping file: QIIME format can have extra meta data columns')
	parser.add_argument('--reads', dest="reads", default='paired', choices=['paired', 'forward'], help='PE or forward reads')
	parser.add_argument('--read_length', type=int, help='Read length, i.e. 2 x 300 bp = 300')
	parser.add_argument('-f','--fwd_primer', dest="F_primer", default='fITS7', help='Forward Primer (fITS7)')
	parser.add_argument('-r','--rev_primer', dest="R_primer", default='ITS4', help='Reverse Primer (ITS4)')
	parser.add_argument('--require_primer', dest="primer", default='on', choices=['on', 'off'], help='Require Fwd primer to be present')
	parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
	parser.add_argument('--barcode_mismatch', default=1, type=int, help='Number of mis-matches allowed in index')
	parser.add_argument('--rescue_forward', default='on', choices=['on', 'off'], help='Rescue Not-merged forward reads')
	parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
	parser.add_argument('--merge_method', default='usearch', choices=['usearch', 'vsearch'], help='Software to use for PE read merging')
	parser.add_argument('-l','--trim_len', default=300, type=int, help='Trim length for reads')
	parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
	parser.add_argument('--full_length', action='store_true', help='Keep only full length reads (no trimming/padding)')
	parser.add_argument('-p','--pad', default='off', choices=['on', 'off'], help='Pad with Ns to a set length')
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH executable')
	parser.add_argument('--sra', action='store_true', help='Input files are from NCBI SRA not direct from illumina')
	parser.add_argument('--cleanup', action='store_true', help='Delete all intermediate files')
	args=parser.parse_args(args)

	#sometimes people add slashes in the output directory, this could be bad, try to fix it
	args.out = re.sub(r'\W+', '', args.out)
			
	#create directory and check for existing logfile
	if not os.path.exists(args.out):
		os.makedirs(args.out)
	
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

	#Now all the data is in folder args.out that needs to be de-multiplexed
	if not args.cpus:
		cpus = multiprocessing.cpu_count()
	else:
		cpus = args.cpus

	#check folder if files are gzipped, then gunzip them
	#try to gunzip files
	gzip_list = []
	for file in os.listdir(args.input):
		if file.endswith(".fastq.gz"):
			gzip_list.append(file)
	if gzip_list:
		amptklib.log.info("Gzipped files detected, uncompressing")
		for file in gzip_list:
			amptklib.log.debug("Uncompressing %s" % file)
			OutName = os.path.join(args.input, os.path.splitext(file)[0])
			amptklib.Funzip(os.path.join(args.input, file), OutName, cpus) 

	#check for mapping file, if exists, then use names from first column only for filenames
	SampleData = {}
	Barcodes = {}
	RevBarcodes = {}
	FwdPrimer = ''
	RevPrimer = ''
	if args.mapping_file:
		if not os.path.isfile(args.mapping_file):
			amptklib.error("Mapping file is not valid: %s" % args.mapping_file)
			sys.exit(1)
		SampleData, Barcodes, RevBarcodes, FwdPrimer, RevPrimer = amptklib.parseMappingFileNEW(args.mapping_file)  
		mapdata = amptklib.parseMappingFileIllumina(args.mapping_file)
		#forward primer in first item in tuple, reverse in second
		sample_names = list(SampleData.keys())
		#loop through the files in the folder and get the ones in the sample_names lit
		filenames = []
		for file in os.listdir(args.input):
			if file.startswith(tuple(sample_names)):
				if file.endswith('.fastq'):
					filenames.append(file)
	
		if len(filenames) < 1:
			amptklib.log.error("Found 0 valid files from mapping file. Mapping file SampleID must match start of filenames")
			sys.exit(1)

	else: #if not then search through and find all the files you can in the folder
		'''get filenames, store in list, Illumina file names look like the following:
		<sample name>_<i5>-<i7>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz'''

		#now get the FASTQ files and proceed
		filenames = []
		for file in os.listdir(args.input):
			if file.endswith(".fastq"):
				filenames.append(file)
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

	#if files are from SRA, then do something different as they are already merged
	if args.sra:
		#take list of filenames, move over to output folder
		sampleDict = {}
		fastq_for = []
		for x in filenames:
			rename = os.path.basename(x).split(".f",-1)[0]
			sampleDict[rename] = 'unknown'
			shutil.copyfile(os.path.join(args.input, x), os.path.join(args.out, rename+'.fq'))
			fastq_for.append(os.path.join(args.out, rename+'.fq'))
		args.reads = 'forward'
	else:
		if len(filenames) % 2 != 0:
			print("Check your input files, they do not seem to be properly paired")
			sys.exit(1)
		#check list for files, i.e. they need to have _R1 and _R2 in the filenames, otherwise throw exception
		if not any('_R1' in x for x in filenames):
			amptklib.log.error("Did not find valid FASTQ files.  Your files must have _R1 and _R2 in filename, rename your files and restart script.")
			sys.exit(1)
		uniq_names = []
		fastq_for = []
		fastq_rev = []
		sampleDict = {}
		map = args.out + '.filenames.txt'
		with open(map, 'w') as map_file:
			map_file.write("Name\t[i5]\t[i7]\tLane\tSet_num\n")
			for item in sorted(filenames):
				if '_R1' in item:
					fastq_for.append(os.path.join(args.input, item))
				if '_R2' in item:
					fastq_rev.append(os.path.join(args.input, item))
				column = item.split("_")
				if column[0] not in uniq_names:
					uniq_names.append(column[0])
					if "-" in column[1]:
						barcode = column[1].split("-")#looking here for the linker between i5 and i7 seqs
						i5 = barcode[0]
						i7 = barcode[1]
						try:
							map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], i5, i7, column[2], column[4].split(".",1)[0]))
						except IndexError:
							amptklib.log.debug("Non-standard names detected, skipping mapping file")              
					else:
						i5 = column[1]
						i7 = "None"
						try:
							map_file.write("%s\t%s\t%s\t%s\t%s\n" % (column[0], i5, i7, column[2], column[4].split(".",1)[0]))
						except IndexError:
							amptklib.log.debug("Non-standard names detected, skipping mapping file")
					if i7 != "None":
						sampleDict[column[0]] = i5+'-'+i7
					else:
						sampleDict[column[0]] = i5

	if args.full_length and args.primer == 'off':
		amptklib.log.info('--full_length is not compatible with --require_primer off, turning --full_length off')
		args.full_length = False

	#tell user about number of cores using
	amptklib.log.info('Demuxing data using {:} cpus'.format(cpus))
	amptklib.log.info('Dropping reads less than {:} bp and setting lossless trimming to {:} bp.'.format(args.min_len, args.trim_len))

	#zip read lists into a single list of tuples
	if args.reads == 'paired':
		amptklib.log.info("Strip Primers and Merge PE reads. FwdPrimer: {:} RevPrimer: {:}".format(FwdPrimer, RevPrimer))
		readList = list(zip(fastq_for, fastq_rev))
		amptklib.runMultiProgress(safe_run, readList, cpus, args=args)
	else:
		amptklib.log.info("Strip Primers. FwdPrimer: {:} RevPrimer: {:}".format(FwdPrimer, RevPrimer))
		amptklib.runMultiProgress(safe_run2, fastq_for, cpus, args=args)
	
	#cleanup to save space
	if gzip_list:
		for file in gzip_list:
			file = file.replace('.gz', '')
			amptklib.removefile(os.path.join(args.input, file))
	print("-------------------------------------------------------")
	#Now concatenate all of the demuxed files together
	amptklib.log.info("Concatenating Demuxed Files")

	catDemux = args.out + '.demux.fq'
	with open(catDemux, 'w') as outfile:
		for filename in glob.glob(os.path.join(args.out,'*.demux.fq')):
			if filename == catDemux:
				continue
			with open(filename, 'r') as readfile:
				shutil.copyfileobj(readfile, outfile)

	#parse the stats
	#(Total, ForPrimerFound, RevPrimerFound, multiHits, TooShort, ValidSeqs))
	finalstats = [0,0,0,0,0,0]
	for file in os.listdir(args.out):
		if file.endswith('.stats'):
			with open(os.path.join(args.out, file), 'r') as statsfile:
				line = statsfile.readline()
				line = line.replace('\n', '')
				newstats = line.split(',')
				newstats = [int(i) for i in newstats]
				for x, num in enumerate(newstats):
					finalstats[x] += num           
	amptklib.log.info('{0:,}'.format(finalstats[0])+' total reads')
	amptklib.log.info('{0:,}'.format(finalstats[1])+' Fwd Primer found, {0:,}'.format(finalstats[2])+ ' Rev Primer found')
	amptklib.log.info('{0:,}'.format(finalstats[3])+ ' discarded Primer incompatibility')
	amptklib.log.info('{0:,}'.format(finalstats[4])+' discarded too short (< %i bp)' % args.min_len)
	amptklib.log.info('{0:,}'.format(finalstats[5])+' valid output reads')


	#now loop through data and find barcoded samples, counting each.....
	BarcodeCount = {}
	with open(catDemux, 'r') as input:
		header = itertools.islice(input, 0, None, 4)
		for line in header:
			ID = line.split("=")[-1].split(";")[0]
			if ID not in BarcodeCount:
				BarcodeCount[ID] = 1
			else:
				BarcodeCount[ID] += 1

	#now let's count the barcodes found and count the number of times they are found.
	barcode_counts = "%30s:  %s" % ('Sample', 'Count')
	for k,v in natsorted(list(BarcodeCount.items()), key=lambda k_v: k_v[1], reverse=True):
		barcode_counts += "\n%30s:  %s" % (k, str(BarcodeCount[k]))
	amptklib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

	genericmapfile = args.out + '.mapping_file.txt'
	if not args.mapping_file:
		#create a generic mappingfile for downstream processes
		amptklib.CreateGenericMappingFileIllumina(sampleDict, FwdPrimer, RevPrimer, genericmapfile, BarcodeCount)
	else:
		amptklib.updateMappingFile(args.mapping_file, BarcodeCount, genericmapfile)
	
	#compress the output to save space
	FinalDemux = catDemux+'.gz'
	amptklib.Fzip(catDemux, FinalDemux, cpus)
	amptklib.removefile(catDemux)

	#get file size
	filesize = os.path.getsize(FinalDemux)
	readablesize = amptklib.convertSize(filesize)
	amptklib.log.info("Output file:  %s (%s)" % (FinalDemux, readablesize))
	amptklib.log.info("Mapping file: %s" % genericmapfile)
	if args.cleanup:
		shutil.rmtree(args.out)
	print("-------------------------------------------------------")
	if 'darwin' in sys.platform:
		print(col.WARN + "\nExample of next cmd: " + col.END + "amptk cluster -i %s -o out\n" % (FinalDemux))
	else:
		print("\nExample of next cmd: amptk cluster -i %s -o out\n" % (FinalDemux))

if __name__ == "__main__":
	main(args)