#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import subprocess
import logging
import shutil
from Bio import SeqIO
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x
    
def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize
    
def main(args):
	parser=argparse.ArgumentParser(prog='amptk-OTU_cluster_ref.py', usage="%(prog)s [options] -i file.demux.fq\n%(prog)s -h for help menu",
		description='''Script runs UPARSE OTU clustering.
		Requires USEARCH by Robert C. Edgar: http://drive5.com/usearch''',
		epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--fastq', dest="FASTQ", required=True, help='FASTQ file (Required)')
	parser.add_argument('-o','--out', help='Base output name')
	parser.add_argument('-e','--maxee', default='1.0', help='Quality trim EE value')
	parser.add_argument('-p','--pct_otu', default='97', help="OTU Clustering Percent")
	parser.add_argument('--id', default='97', help="Threshold for alignment")
	parser.add_argument('-m','--minsize', default='2', help='Min identical seqs to process')
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
	parser.add_argument('--map_filtered', action='store_true', help='map quality filtered reads back to OTUs')
	parser.add_argument('-d','--db', required=True, help='Reference Database [ITS,ITS1,ITS2,16S,LSU,COI,custom]')
	parser.add_argument('--utax_db', help='UTAX Reference Database')
	parser.add_argument('--utax_cutoff', default=0.8, type=restricted_float, help='UTAX confidence value threshold.')
	parser.add_argument('--utax_level', default='k', choices=['k','p','c','o','f','g','s'], help='UTAX classification level to retain')
	parser.add_argument('--mock', default='synmock', help='Spike-in mock community (fasta)')
	parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
	parser.add_argument('--closed_ref_only', action='store_true', help='Only run closed reference clustering')
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


	taxonomyLookup = {'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class', 'o': 'Order', 'f': 'Family', 'g': 'Genus', 's': 'Species'}

	#remove logfile if exists
	log_name = base + '.amptk-cluster_ref.log'
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

	#Setup DB locations and names, etc
	DBdir = os.path.join(parentdir, 'DB')   
	DataBase = {'ITS1': (os.path.join(DBdir, 'ITS.udb'), os.path.join(DBdir, 'ITS1_UTAX.udb')), 
				'ITS2': (os.path.join(DBdir, 'ITS.udb'), os.path.join(DBdir, 'ITS2_UTAX.udb')), 
				'ITS': (os.path.join(DBdir, 'ITS.udb'), os.path.join(DBdir, 'ITS_UTAX.udb')), 
				'16S': (os.path.join(DBdir, '16S.udb'), os.path.join(DBdir, '16S.udb')), 
				'LSU': (os.path.join(DBdir, 'LSU.udb'), os.path.join(DBdir, 'LSU_UTAX.udb')), 
				'COI': (os.path.join(DBdir, 'COI.udb'), os.path.join(DBdir, 'COI_UTAX.udb'))}

	#setup refDB
	amptklib.log.info("Checking Reference Database")
	if args.db in DataBase:
		#need to write to fasta from vsearch UDB
		DB = os.path.join(tmp, args.db+'.extracted.fa')
		cmd = ['vsearch', '--udb2fasta', DataBase.get(args.db)[0], '--output', DB]
		amptklib.runSubprocess(cmd, amptklib.log)
	else:
		DB = os.path.abspath(args.db)
	refDB = os.path.join(tmp, 'reference_DB.fa')
	if args.mock:
		if args.mock == 'synmock':
			mock = os.path.join(parentdir, 'DB', 'amptk_synmock.fa')
		else:
			mock = os.path.abspath(args.mock)
	seen = []
	with open(refDB, 'w') as output:
		if args.mock:
			with open(mock) as input1:
				for rec in SeqIO.parse(input1, 'fasta'):
					if not rec.id in seen:
						SeqIO.write(rec, output, 'fasta')
					else:
						amptklib.log.error("Duplicate ID's in Ref DB: %s, exiting" % rec.id)
						sys.exit(1)
		with open(DB) as input2:
			for rec in SeqIO.parse(input2, 'fasta'):
				if not rec.id in seen:
					SeqIO.write(rec, output, 'fasta')
				else:
					amptklib.log.error("Duplicate ID's in Ref DB: %s, exiting" % rec.id)
					sys.exit(1)

	#get utax_database
	if args.db in DataBase:
		utaxDB = DataBase.get(args.db)[1]
	else:
		if not args.closed_ref_only:
			if args.utax_db:
				utaxDB = os.path.abspath(args.utax_db)
			else:
				amptklib.log.error("%s not pre-installed DB, must then also specify valid UTAX database via --utax_db" % args.db)
				sys.exit(1)

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
	qtrimtotal = amptklib.countfastq(filter_out)
	amptklib.log.info('{0:,}'.format(qtrimtotal) + ' reads passed')
	
	#now run full length dereplication
	derep_out = os.path.join(tmp, base + '.EE' + args.maxee + '.derep.fa')
	amptklib.log.info("De-replication (remove duplicate reads)")
	cmd = ['vsearch', '--derep_fulllength', filter_fasta, '--sizeout', '--output', derep_out, '--threads', str(cpus), '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfasta(derep_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads passed')

	#now run sort by size
	sort_out = os.path.join(tmp, base + '.EE' + args.maxee + '.sort.fa')
	amptklib.log.info("Sorting reads by size: removing reads seen less than %s times" % args.minsize)
	cmd = ['vsearch', '--sortbysize', derep_out, '--minsize', args.minsize, '--output', sort_out, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfasta(sort_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads passed')

	#chimera detection
	#first run through de novo chimera detection
	amptklib.log.info("De novo chimera detection (VSEARCH)")
	chimera_out = os.path.join(tmp, base + '.EE' + args.maxee + '.chimera_check.fa')
	cmd = ['vsearch', '--uchime_denovo', sort_out, '--relabel', 'Seq', '--sizeout', '--nonchimeras', chimera_out, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfasta(chimera_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads passed')
  
	#now run uchime_ref
	uchime_out = os.path.join(tmp, base + '.EE' + args.maxee + '.uchime.otus.fa')
	#now run chimera filtering if all checks out
	amptklib.log.info("Chimera Filtering (VSEARCH)")
	cmd = ['vsearch', '--mindiv', '1.0', '--uchime_ref', chimera_out, '--db', refDB, '--sizeout', '--nonchimeras', uchime_out, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)
	total = amptklib.countfasta(uchime_out)
	amptklib.log.info('{0:,}'.format(total) + ' OTUs passed')
	
	#now run usearch_global versus reference database
	align_out = os.path.join(tmp, base + '.align.uc')
	pident = int(args.id) * 0.01
	amptklib.log.info("Reference Clustering using Global Alignment, %s%% identity" % args.id)
	cmd = ['vsearch', '--usearch_global', uchime_out, '--db', refDB, '--id', str(pident), '--output_no_hits', '--top_hits_only', '--notrunclabels', '--uc', align_out, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)

	#parse results
	ref_results = {}
	nohits = []
	with open(align_out, 'r') as alignment:
		for line in alignment:
			line = line.replace('\n', '')
			col = line.split('\t')
			counts = col[8].split(';')
			counts = int(counts[1].replace('size=', ''))
			if col[3] == '*':
				nohits.append(col[8])
				continue
			if float(col[3]) >= float(args.id):
				if not col[8] in ref_results:
					ref_results[col[8]] = (col[9], col[3], counts)
				else:
					print("Error: %s duplicated ID" % col[8])
			else:
				nohits.append(col[8])

	#summarize results from first ref clustering
	num_refcluster = len(ref_results)
	seqs_refcluster = 0
	for k,v in list(ref_results.items()):
		seqs_refcluster += v[2]
	amptklib.log.info("%i OTUs classified " % num_refcluster + "({0:.0f}%".format(seqs_refcluster/float(qtrimtotal)* 100)+ " of reads)")

	#get ref clustered hits to file with taxonomy
	ref_clustered = os.path.join(tmp, base+'.ref_clustered.fa')
	with open(ref_clustered, 'w') as refoutput:
		with open(uchime_out, 'r') as input:
			otu_counter = 1
			for rec in SeqIO.parse(input, 'fasta'):
				if rec.id in ref_results:
					res = ref_results.get(rec.id)
					pident = res[1]
					tax = res[0]
					newID = 'OTU'+str(otu_counter)+';pident='+pident+';'+tax
					rec.id = newID
					rec.name = ''
					rec.description = ''
					SeqIO.write(rec, refoutput, 'fasta')
					otu_counter += 1

	if not args.closed_ref_only:
		#get nohits file to run clustering
		utax_ref = os.path.join(tmp, base + '.EE' + args.maxee + '.utax_ref.fa')
		with open(utax_ref, 'w') as output:
			with open(uchime_out, 'r') as input:
				for rec in SeqIO.parse(input, 'fasta'):
					if rec.id in nohits:
						SeqIO.write(rec, output, 'fasta')

		#input needs to be sorted, so 
		ref_sort = os.path.join(tmp, base+'.utax_ref.sorted.fa')
		cmd = ['vsearch', '--sortbysize', utax_ref, '--minsize', args.minsize, '--output', ref_sort, '--threads', str(cpus)]
		amptklib.runSubprocess(cmd, amptklib.log)
		   
		#now run clustering algorithm on those not found in reference database
		radius = str(100 - int(args.pct_otu))
		otu_out = os.path.join(tmp, base + '.EE' + args.maxee + '.otus.fa')
		amptklib.log.info("De novo Clustering remaining sequences (UPARSE)")
		cmd = [usearch, '-cluster_otus', ref_sort, '-relabel', 'OTU', '-otu_radius_pct', radius, '-otus', otu_out]
		amptklib.runSubprocess(cmd, amptklib.log)
		total = amptklib.countfasta(otu_out)
		amptklib.log.info('{0:,}'.format(total) + ' de novo OTUs')

		#try utax reference clustering
		amptklib.log.info("Reference Clustering de novo OTUs using UTAX")
		cmd = [usearch, '-cluster_otus_utax', otu_out, '-db', utaxDB, '-utax_cutoff', str(args.utax_cutoff), '-utax_level', 's', '-strand', 'plus', '-utaxout', os.path.join(tmp, base+'.utax.out')]
		amptklib.runSubprocess(cmd, amptklib.log)
		#setup tax filtering
		tax_values = ['k','p','c','o','f','g','s']
		filter_index = tax_values.index(args.utax_level)
		filt_tax_values = [s + ':' for s in tax_values[filter_index:]]
		#get results from utax
		with open(ref_clustered, 'a') as output:
			seqDict = SeqIO.index(otu_out, 'fasta')
			utaxresults = []
			with open(os.path.join(tmp, base+'.utax.out'), 'r') as utax:
				for line in utax:
					line = line.replace('\n', '')
					col = line.split('\t')
					ID = col[0]
					tax = col[2]
					if any(x in tax for x in filt_tax_values):
						record = seqDict[ID]
						record.id = 'OTU'+str(otu_counter)+';UTAX;tax='+tax
						record.name = ''
						record.description = ''
						SeqIO.write(record, output, 'fasta')
						otu_counter += 1
		total = amptklib.countfasta(ref_clustered) - num_refcluster
		amptklib.log.info('{0:,}'.format(total) + ' classified to %s' % taxonomyLookup.get(args.utax_level))

	#clean up padded N's
	amptklib.log.info("Cleaning up padding from OTUs")
	otu_clean = os.path.join(tmp, base + '.clean.otus.fa')
	amptklib.fasta_strip_padding(ref_clustered, otu_clean)           
	total = amptklib.countfasta(otu_clean)
	amptklib.log.info('{0:,}'.format(total) + ' total OTUs')
	   
	#now map reads back to OTUs
	uc_out = os.path.join(tmp, base + '.EE' + args.maxee + '.mapping.uc')
	otu_table = os.path.join(tmp, base + '.EE' + args.maxee + '.otu_table.txt')
	#setup reads to map
	if args.map_filtered:
		reads = filter_fasta
	else:
		reads = orig_fasta
	amptklib.log.info("Mapping Reads to OTUs and Building OTU table")
	cmd = ['vsearch', '--usearch_global', reads, '--strand', 'plus', '--id', '0.97', '--db', otu_clean, '--uc', uc_out, '--otutabout', otu_table, '--threads', str(cpus)]
	amptklib.runSubprocess(cmd, amptklib.log)

	#count reads mapped
	total = amptklib.line_count2(uc_out)
	amptklib.log.info('{0:,}'.format(total) + ' reads mapped to OTUs '+ '({0:.0f}%)'.format(total/float(orig_total)* 100))

	#Move files around, delete tmp if argument passed.
	currentdir = os.getcwd()
	final_otu = os.path.join(currentdir, base + '.cluster.otus.fa')
	shutil.copyfile(otu_clean, final_otu)
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