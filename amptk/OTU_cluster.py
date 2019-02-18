#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import inspect
import shutil
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-OTU_cluster.py', usage="%(prog)s [options] -i file.demux.fq\n%(prog)s -h for help menu",
		description='''Script runs UPARSE OTU clustering.
		Requires USEARCH by Robert C. Edgar: http://drive5.com/usearch''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
	parser.add_argument('-o','--out', help='Base output name')
	parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
	parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
	parser.add_argument('-m','--minsize', default='2', help='Min size to keep for clustering')
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
	parser.add_argument('--uchime_ref', help='Run UCHIME REF [ITS,16S,LSU,COI,custom]')
	parser.add_argument('--map_filtered', action='store_true', help='map quality filtered reads back to OTUs')
	parser.add_argument('--unoise', action='store_true', help='Run De-noising (UNOISE)')
	parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
	parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
	args=parser.parse_args(args)
	
	parentdir = os.path.join(os.path.dirname(amptklib.__file__))

	#get basename if not args.out passed
	if args.out:
		base = args.out
	else:
		if 'demux' in args.FASTQ:
			base = os.path.basename(args.FASTQ).split('.demux')[0]
		else:
			base = os.path.basename(args.FASTQ).split('.f')[0]

	#remove logfile if exists
	log_name = base + '.amptk-cluster.log'
	if os.path.isfile(log_name):
		os.remove(log_name)

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

	#get number of cpus
	if args.cpus:
		cpus = args.cpus
	else:
		cpus = amptklib.getCPUS()

	#make tmp folder
	tmp = base + '_tmp'
	if not os.path.exists(tmp):
		os.makedirs(tmp)

	#Count FASTQ records
	amptklib.log.info("Loading FASTQ Records")
	#convert to FASTA for mapping
	orig_fasta = os.path.join(tmp, base+'.orig.fa')
	cmd = ['vsearch', '--fastq_filter', args.FASTQ, '--fastaout', orig_fasta, '--fastq_qmax', '55', '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	orig_total = amptklib.countfasta(orig_fasta)
	size = amptklib.checkfastqsize(args.FASTQ)
	readablesize = amptklib.convertSize(size)
	amptklib.log.info('{0:,}'.format(orig_total) + ' reads (' + readablesize + ')')

	#Expected Errors filtering step
	filter_out = os.path.join(tmp, base + '.EE' + args.maxee + '.filter.fq')
	filter_fasta = os.path.join(tmp, base + '.EE' + args.maxee + '.filter.fa')
	amptklib.log.info("Quality Filtering, expected errors < %s" % args.maxee)
	cmd = ['vsearch', '--fastq_filter', args.FASTQ, '--fastq_maxee', str(args.maxee), '--fastqout', filter_out, '--fastaout', filter_fasta, '--fastq_qmax', '55', '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfastq(filter_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads passed')

	#now run full length dereplication
	derep_out = os.path.join(tmp, base + '.EE' + args.maxee + '.derep.fa')
	amptklib.log.info("De-replication (remove duplicate reads)")
	cmd = ['vsearch', '--derep_fulllength', filter_fasta, '--sizeout', '--output', derep_out, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfasta(derep_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads passed')

	#optional run UNOISE
	if args.unoise:
		unoise_out = unoise_out = os.path.join(tmp, base + '.EE' + args.maxee + '.denoised.fa')
		amptklib.log.info("Denoising Data with UNOISE")
		cmd = [usearch, '-cluster_fast', derep_out, '-centroids', unoise_out, '-id', '0.9', '--maxdiffs', '5', '-abskew', '10', '-sizein', '-sizeout', '-sort', 'size', '-threads', str(cpus)]
		amptklib.runSubprocess(cmd, amptklib.log)   
		total = amptklib.countfasta(unoise_out)
		amptklib.log.info('{0:,}'.format(total) + ' reads passed')
	else:
		unoise_out = derep_out

	#now sort by size remove singletons
	sort_out = os.path.join(tmp, base + '.EE' + args.maxee + '.sort.fa')
	cmd = ['vsearch', '--sortbysize', unoise_out, '--minsize', args.minsize, '--output', sort_out, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)

	#now run clustering algorithm
	radius = str(100 - int(args.pct_otu))
	otu_out = os.path.join(tmp, base + '.EE' + args.maxee + '.otus.fa')
	amptklib.log.info("Clustering OTUs (UPARSE)") 
	cmd = [usearch, '-cluster_otus', sort_out, '-relabel', 'OTU', '-otu_radius_pct', radius, '-otus', otu_out, '-threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	numOTUs = amptklib.countfasta(otu_out)
	amptklib.log.info('{0:,}'.format(numOTUs) + ' OTUs')

	#clean up padded N's
	amptklib.log.info("Cleaning up padding from OTUs")
	otu_clean = os.path.join(tmp, base + '.EE' + args.maxee + '.clean.otus.fa')
	amptklib.fasta_strip_padding(otu_out, otu_clean)

	#optional UCHIME Ref
	if not args.uchime_ref:
		uchime_out = otu_clean
	else:
		uchime_out = os.path.join(tmp, base + '.EE' + args.maxee + '.uchime.otus.fa')
		#check if file is present, remove from previous run if it is.
		if os.path.isfile(uchime_out):
			os.remove(uchime_out)
		#R. Edgar now says using largest DB is better for UCHIME, so use the one distributed with taxonomy
		if args.uchime_ref in ['ITS', '16S', 'LSU', 'COI']: #test if it is one that is setup, otherwise default to full path
			uchime_db = os.path.join(parentdir, 'DB', args.uchime_ref+'.udb')
			if not os.path.isfile(uchime_db):
				amptklib.log.error("Database not properly configured, run `amptk install` to setup DB, skipping chimera filtering")
				uchime_out = otu_clean
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
				uchime_out = otu_clean
		#now run chimera filtering if all checks out
		if not os.path.isfile(uchime_out):
			amptklib.log.info("Chimera Filtering (VSEARCH) using %s DB" % args.uchime_ref)
			cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', otu_clean, '--db', uchime_db, '--nonchimeras', uchime_out, '--threads', str(cpus)]
			amptklib.runSubprocess(cmd, amptklib.log)
			total = amptklib.countfasta(uchime_out)
			uchime_chimeras = numOTUs - total
			amptklib.log.info('{0:,}'.format(total) + ' OTUs passed, ' + '{0:,}'.format(uchime_chimeras) + ' ref chimeras')

	#Filter out OTUs in wrong orientation
	amptklib.log.info('Validating OTU orientation')
	passingOTUs = os.path.join(tmp, base+'.passed.otus.fa')
	numKept, numDropped = amptklib.validateorientation(tmp, sort_out, uchime_out, passingOTUs)
	amptklib.log.info('{:,} OTUs validated ({:,} dropped)'.format(numKept, numDropped))

	#now map reads back to OTUs and build OTU table
	uc_out = os.path.join(tmp, base + '.EE' + args.maxee + '.mapping.uc')
	otu_table = os.path.join(tmp, base + '.EE' + args.maxee + '.otu_table.txt')
	#setup reads to map
	if args.map_filtered:
		reads = filter_fasta
	else:
		reads = orig_fasta
	amptklib.log.info("Mapping Reads to OTUs and Building OTU table")
	cmd = ['vsearch', '--usearch_global', reads, '--strand', 'plus', '--id', '0.97', '--db', passingOTUs, '--uc', uc_out, '--otutabout', otu_table, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)

	#count reads mapped
	total = amptklib.line_count2(uc_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

	#Move files around, delete tmp if argument passed.
	currentdir = os.getcwd()
	final_otu = os.path.join(currentdir, base + '.cluster.otus.fa')
	shutil.copyfile(passingOTUs, final_otu)
	final_otu_table = os.path.join(currentdir, base + '.otu_table.txt')
	shutil.copyfile(otu_table, final_otu_table)
	if not args.debug:
		shutil.rmtree(tmp)

	#Print location of files to STDOUT
	print("-------------------------------------------------------")
	print("OTU Clustering Script has Finished Successfully")
	print("-------------------------------------------------------")
	if not not args.debug:
		print("Tmp Folder of files: %s" % tmp)
	print("Clustered OTUs: %s" % os.path.basename(final_otu))
	print("OTU Table: %s" % os.path.basename(final_otu_table))
	print("-------------------------------------------------------")

	otu_print = final_otu.split('/')[-1]
	tab_print = final_otu_table.split('/')[-1]
	if 'darwin' in sys.platform:
		print(colr.WARN + "\nExample of next cmd:" + colr.END + " amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print))
	else:
		print("\nExample of next cmd: amptk filter -i %s -f %s -b <mock barcode>\n" % (tab_print, otu_print))

if __name__ == "__main__":
	main(args)
