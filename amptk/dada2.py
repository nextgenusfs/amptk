#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import sys
import os
import argparse
import logging
import shutil
import subprocess
import numpy as np
from natsort import natsorted
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'
    
def folder2list(input, ending):
    names = []
    if not os.path.isdir(input):
        return False
    else:
        for x in os.listdir(input):
            if x.endswith(ending):
                x = os.path.join(input, x)
                names.append(x)
    return names

def splitDemux2(input, outputdir, args=False):
    for title, seq, qual in FastqGeneralIterator(open(input)):
        sample = title.split('barcodelabel=')[1].split(';')[0]
        sample = sample.replace(';', '')
        if not args.length:
            with open(os.path.join(outputdir, sample+'.fastq'), 'a') as output:
                output.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        else:
            if len(seq) >= int(args.length):
                with open(os.path.join(outputdir, sample+'.fastq'), 'a') as output:
                    output.write("@%s\n%s\n+\n%s\n" % (title, seq[:int(args.length):], qual[:int(args.length)]))

def getAvgLength(input):
    AvgLength = []
    for title, seq, qual in FastqGeneralIterator(open(input)):
        AvgLength.append(len(seq))
    Average = sum(AvgLength) / float(len(AvgLength))
    Min = min(AvgLength)
    Max = max(AvgLength)
    a = np.array(AvgLength)
    nintyfive = np.percentile(a, 5)
    return (Average, Min, Max, int(nintyfive))

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-dada2.py',
		description='''Script takes output from amptk pre-processing and runs DADA2''',
		epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--fastq', required=True, help='Input Demuxed containing FASTQ')
	parser.add_argument('-o','--out', help='Output Basename')
	parser.add_argument('-m','--min_reads', default=10, type=int, help="Minimum number of reads after Q filtering to run DADA2 on")
	parser.add_argument('-l','--length', type=int, help='Length to truncate reads')
	parser.add_argument('-e','--maxee', default='1.0', help='MaxEE quality filtering')
	parser.add_argument('-p','--pct_otu', default='97', help="Biological OTU Clustering Percent")
	parser.add_argument('--platform', default='ion', choices=['ion', 'illumina', '454'], help='Sequencing platform')
	parser.add_argument('--chimera_method', default='consensus', choices=['consensus', 'pooled', 'per-sample'], help='bimera removal method')
	parser.add_argument('--uchime_ref', help='Run UCHIME REF [ITS,16S,LSU,COI,custom]')
	parser.add_argument('--pool', action='store_true', help='Pool all sequences together for DADA2')
	parser.add_argument('--debug', action='store_true', help='Keep all intermediate files')
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
	parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
	args=parser.parse_args(args)

	parentdir = os.path.join(os.path.dirname(amptklib.__file__))
	dada2script = os.path.join(parentdir, 'dada2_pipeline_nofilt.R')

	#get basename if not args.out passed
	if args.out:
		base = args.out
	else:
		if 'demux' in args.fastq:
			base = os.path.basename(args.fastq).split('.demux')[0]
		else:
			base = os.path.basename(args.fastq).split('.f')[0]

	#remove logfile if exists
	log_name = base + '.amptk-dada2.log'
	if os.path.isfile(log_name):
		amptklib.removefile(log_name)

	amptklib.setupLogging(log_name)
	FNULL = open(os.devnull, 'w')
	cmd_args = " ".join(sys.argv)+'\n'
	amptklib.log.debug(cmd_args)
	print("-------------------------------------------------------")
	#initialize script, log system info and usearch version
	amptklib.SystemInfo()
	#Do a version check
	usearch = args.usearch
	amptklib.versionDependencyChecks(usearch)

	#get number of cores
	if args.cpus:
		CORES = str(args.cpus)
	else:
		CORES = str(amptklib.getCPUS())

	#check dependencies
	programs = ['Rscript']
	amptklib.CheckDependencies(programs)
	Rversions = amptklib.checkRversion()
	R_pass = '3.2.1'
	dada2_pass = '1.3.3'

	#check dada2 first, if good move on, otherwise issue warning
	if not amptklib.gvc(Rversions[1], dada2_pass):
		amptklib.log.error("R v%s; DADA2 v%s detected, need atleast v%s" % (Rversions[0], Rversions[1], dada2_pass))
		amptklib.log.error("See: http://benjjneb.github.io/dada2/dada-installation.html")
		sys.exit(1)
	amptklib.log.info("R v%s; DADA2 v%s" % (Rversions[0], Rversions[1]))

	#Count FASTQ records and remove 3' N's as dada2 can't handle them
	amptklib.log.info("Loading FASTQ Records")
	no_ns = base+'.cleaned_input.fq'
	if args.fastq.endswith('.gz'):
		fastqInput = args.fastq.replace('.gz', '')
		amptklib.Funzip(os.path.abspath(args.fastq), os.path.basename(fastqInput), CORES)
	else:
		fastqInput = os.path.abspath(args.fastq)
	amptklib.fastq_strip_padding(os.path.basename(fastqInput), no_ns)
	demuxtmp = base+'.original.fa'
	cmd = ['vsearch', '--fastq_filter', os.path.abspath(no_ns),'--fastq_qmax', '55', '--fastaout', demuxtmp, '--threads', CORES]
	amptklib.runSubprocess(cmd, amptklib.log)
	orig_total = amptklib.countfasta(demuxtmp)
	size = amptklib.checkfastqsize(no_ns)
	readablesize = amptklib.convertSize(size)
	amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

	#quality filter
	amptklib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
	derep = base+'.qual-filtered.fq'
	filtercmd = ['vsearch', '--fastq_filter', no_ns, '--fastq_maxee', str(args.maxee), '--fastqout', derep, '--fastq_qmax', '55', '--fastq_maxns', '0', '--threads', CORES]
	amptklib.runSubprocess(filtercmd, amptklib.log)
	total = amptklib.countfastq(derep)
	amptklib.log.info('{0:,}'.format(total) + ' reads passed')

	#split into individual files
	amptklib.log.info("Splitting FASTQ file by Sample into individual files")
	filtfolder = base+'_filtered'
	if os.path.isdir(filtfolder):
		shutil.rmtree(filtfolder)
	os.makedirs(filtfolder)
	splitDemux2(derep, filtfolder, args=args)

	#check for minimum number of reads in each sample
	remove = []
	files = [i for i in os.listdir(filtfolder) if i.endswith('.fastq')]
	for x in files:
		if amptklib.countfastq(os.path.join(filtfolder, x)) < args.min_reads:
			remove.append(x)
	if len(remove) > 0:
		amptklib.log.info("Dropping %s as fewer than %i reads" % (', '.join(remove), args.min_reads))
		for y in remove:
			os.remove(os.path.join(filtfolder, y))

	#now run DADA2 on filtered folder
	amptklib.log.info("Running DADA2 pipeline")
	dada2log = base+'.dada2.Rscript.log'
	dada2out = base+'.dada2.csv'
	#check pooling vs notpooled, default is not pooled.
	if args.pool:
		POOL = 'TRUE'
	else:
		POOL = 'FALSE'
	with open(dada2log, 'w') as logfile:
		subprocess.call(['Rscript', '--vanilla', dada2script, filtfolder, dada2out, args.platform, POOL, CORES, args.chimera_method], stdout = logfile, stderr = logfile)

	#check for results
	if not os.path.isfile(dada2out):
		amptklib.log.error("DADA2 run failed, please check %s logfile" % dada2log)
		sys.exit(1)
	
	#now process the output, pull out fasta, rename, etc
	fastaout = base+'.otus.tmp'
	OTUCounts = {}
	counter = 1
	with open(fastaout, 'w') as writefasta:
		with open(dada2out, 'r') as input:
			next(input)
			for line in input:
				line = line.replace('\n', '')
				line = line.replace('"', '')
				cols = line.split(',')
				Seq = cols[0]
				countList = [int(x) for x in cols[1:]]
				counts = sum(countList)
				ID = 'ASV'+str(counter)
				if not ID in OTUCounts:
					OTUCounts[ID] = counts
				writefasta.write(">%s\n%s\n" % (ID, Seq))
				counter += 1

	#get number of bimeras from logfile
	with open(dada2log, 'r') as bimeracheck:
		for line in bimeracheck:
			if line.startswith('Identified '):
				bimeraline = line.split(' ')
				bimeras = int(bimeraline[1])
				totalSeqs = int(bimeraline[5])
	validSeqs = totalSeqs - bimeras
	amptklib.log.info('{0:,}'.format(totalSeqs) + ' total amplicon sequence variants (ASVs)')
	amptklib.log.info('{0:,}'.format(bimeras) + ' denovo chimeras removed')
	amptklib.log.info('{0:,}'.format(validSeqs) + ' valid ASVs')

	#optional UCHIME Ref
	uchime_out = base+'.nonchimeras.fa'
	chimeraFreeTable = base+'.otu_table.txt'
	iSeqs = base+'.ASVs.fa'
	if not args.uchime_ref:
		os.rename(fastaout, iSeqs)
	else:
		#check if file is present, remove from previous run if it is.
		if os.path.isfile(iSeqs):
			amptklib.removefile(iSeqs)
		#R. Edgar now says using largest DB is better for UCHIME, so use the one distributed with taxonomy
		if args.uchime_ref in ['ITS', '16S', 'LSU', 'COI']: #test if it is one that is setup, otherwise default to full path
			uchime_db = os.path.join(parentdir, 'DB', args.uchime_ref+'.udb')
			if not os.path.isfile(uchime_db):
				amptklib.log.error("Database not properly configured, run `amptk install` to setup DB, skipping chimera filtering")
				uchime_out = fastaout
			#since uchime cannot work with udb database, need to extract fasta sequences, do this if 
			if not amptklib.checkfile(os.path.join(parentdir, 'DB',args.uchime_ref+'.extracted.fa')):
				uchime_db = os.path.join(parentdir, 'DB',args.uchime_ref+'.extracted.fa')
				cmd = ['vsearch', '--udb2fasta', os.path.join(parentdir, 'DB', args.uchime_ref+'.udb'), '--output', uchime_db]
				amptklib.runSubprocess(cmd, amptklib.log)
			else:
				uchime_db = os.path.join(parentdir, 'DB',args.uchime_ref+'.extracted.fa')
		else:
			if os.path.isfile(args.uchime_ref):
				uchime_db = os.path.abspath(args.uchime_ref)
			else:
				amptklib.log.error("%s is not a valid file, skipping reference chimera filtering" % args.uchime_ref)
				iSeqs = fastaout
		#now run chimera filtering if all checks out
		if not os.path.isfile(iSeqs):
			amptklib.log.info("Chimera Filtering (VSEARCH) using %s DB" % args.uchime_ref)
			cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', fastaout, '--db', uchime_db, '--nonchimeras', iSeqs, '--threads', CORES]
			amptklib.runSubprocess(cmd, amptklib.log)
			total = amptklib.countfasta(iSeqs)
			uchime_chimeras = validSeqs - total
			amptklib.log.info('{0:,}'.format(total) + ' ASVs passed, ' + '{0:,}'.format(uchime_chimeras) + ' ref chimeras removed')
			if os.path.isfile(fastaout):
				amptklib.removefile(fastaout)

	#setup output files
	dadademux = base+'.dada2.map.uc'
	bioSeqs = base+'.cluster.otus.fa'
	bioTable = base+'.cluster.otu_table.txt'
	uctmp = base+'.map.uc'
	ClusterComp = base+'.ASVs2clusters.txt'

	#Filter out ASVs in wrong orientation
	amptklib.log.info('Validating ASV orientation')
	os.rename(iSeqs, iSeqs+'.bak')
	numKept, numDropped = amptklib.validateorientationDADA2(OTUCounts, iSeqs+'.bak', iSeqs)
	amptklib.log.info('{:,} ASVs validated ({:,} dropped)'.format(numKept, numDropped))
	amptklib.SafeRemove(iSeqs+'.bak')

	#map reads to DADA2 OTUs
	amptklib.log.info("Mapping reads to DADA2 ASVs")
	cmd = ['vsearch', '--usearch_global', demuxtmp, '--db', iSeqs, '--id', '0.97', '--uc', dadademux, '--strand', 'plus', '--otutabout', chimeraFreeTable, '--threads', CORES]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.line_count2(dadademux)
	amptklib.log.info('{0:,}'.format(total) + ' reads mapped to ASVs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

	#cluster
	amptklib.log.info("Clustering ASVs at %s%% to generate biological OTUs" % args.pct_otu)
	radius = float(args.pct_otu) / 100.
	cmd = ['vsearch', '--cluster_smallmem', iSeqs, '--centroids', bioSeqs, '--id', str(radius), '--strand', 'plus', '--relabel', 'OTU', '--qmask', 'none', '--usersort', '--threads', CORES]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfasta(bioSeqs)
	amptklib.log.info('{0:,}'.format(total) + ' OTUs generated')

	#determine where iSeqs clustered
	iSeqmap = base+'.ASV_map.uc'
	cmd = ['vsearch', '--usearch_global', iSeqs, '--db', bioSeqs, '--id', str(radius), '--uc', iSeqmap, '--strand', 'plus', '--threads', CORES]
	amptklib.runSubprocess(cmd, amptklib.log)
	iSeqMapped = {}
	with open(iSeqmap, 'r') as mapping:
		for line in mapping:
			line = line.replace('\n', '')
			cols = line.split('\t')
			OTU = cols[9]
			Hit = cols[8]
			if not OTU in iSeqMapped:
				iSeqMapped[OTU] = [Hit]
			else:
				iSeqMapped[OTU].append(Hit)
	with open(ClusterComp, 'w') as clusters:
		clusters.write('OTU\tASVs\n')
		for k,v in natsorted(list(iSeqMapped.items())):
			clusters.write('%s\t%s\n' % (k, ', '.join(v)))
	#create OTU table
	amptklib.log.info("Mapping reads to OTUs")
	cmd = ['vsearch', '--usearch_global', demuxtmp, '--db', bioSeqs, '--id', '0.97', '--uc', uctmp, '--strand', 'plus', '--otutabout', bioTable, '--threads', CORES]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.line_count2(uctmp)
	amptklib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

	if not args.debug:
		amptklib.removefile(no_ns)
		shutil.rmtree(filtfolder)
		amptklib.removefile(dada2out)
		amptklib.removefile(derep)
		amptklib.removefile(demuxtmp)
		amptklib.removefile(uctmp)
		amptklib.removefile(iSeqmap)
		amptklib.removefile(dadademux)

	#Print location of files to STDOUT
	print("-------------------------------------------------------")
	print("DADA2 Script has Finished Successfully")
	print("-------------------------------------------------------")
	if args.debug:
		print("Tmp Folder of files: %s" % filtfolder)
	print("Amplicon sequence variants: %s" % iSeqs)
	print("ASV OTU Table: %s" % chimeraFreeTable)
	print("Clustered OTUs: %s" % bioSeqs)
	print("OTU Table: %s" % bioTable)
	print("ASVs 2 OTUs: %s" % ClusterComp)
	print("-------------------------------------------------------")

	otu_print = bioSeqs.split('/')[-1]
	tab_print = bioTable.split('/')[-1]
	if 'darwin' in sys.platform:
		print(colr.WARN + "\nExample of next cmd:" + colr.END + " amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print))
	else:
		print("\nExample of next cmd: amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print))   

if __name__ == "__main__":
	main(args)
		