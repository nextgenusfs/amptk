#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import logging
import shutil
import subprocess
from Bio import SeqIO
from amptk import amptklib


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-lulu.py',
		description='''Script runs OTU table post processing LULU to identify low abundance error OTUs''',
		epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--otu_table', required=True, help='Input OTU table')
	parser.add_argument('-f','--fasta', required=True, help='Input OTUs (multi-fasta)')
	parser.add_argument('-o','--out', help='Output folder basename')
	parser.add_argument('--min_ratio_type', default='min', choices=['min', 'avg'], help="LULU minimum ratio threshold")
	parser.add_argument('--min_ratio', default=1, type=int, help="LULU minimum ratio")
	parser.add_argument('--min_match', default=84, type=int, help="LULU minimum match percent identity")
	parser.add_argument('--min_relative_cooccurence', default=95, type=int, help="LULU minimum relative cooccurance")
	parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
	args=parser.parse_args(args)

	#get location of R script
	parentdir = os.path.join(os.path.dirname(amptklib.__file__))
	luluScript = os.path.join(parentdir, 'runLULU.R')

	if not args.out:
		#get base name of files
		if 'otu_table' in args.otu_table:
			base = os.path.basename(args.otu_table).split(".otu_table")[0]
		elif 'final.txt' in args.otu_table:
			base = os.path.basename(args.otu_table).split(".final")[0]
		else:
			base = os.path.basename(args.otu_table).split(".txt")[0]
	else:
		base = args.out

	#remove logfile if exists
	log_name = base + '.amptk-lulu.log'
	if os.path.isfile(log_name):
		amptklib.removefile(log_name)

	amptklib.setupLogging(log_name)
	FNULL = open(os.devnull, 'w')
	cmd_args = " ".join(sys.argv)+'\n'
	amptklib.log.debug(cmd_args)
	print("-------------------------------------------------------")
	#initialize script, log system info and usearch version
	amptklib.SystemInfo()
	amptklib.versionDependencyChecks('usearch9')
	#check dependencies
	programs = ['Rscript', 'vsearch']
	amptklib.CheckDependencies(programs)
	Rversions = amptklib.checkRversion()
	if Rversions[3] == '0.0.0':
		amptklib.log.info("R v%s installed, LULU not installed")
		sys.exit(1)
	else:
		amptklib.log.info("R v%s; LULU v%s" % (Rversions[0], Rversions[3]))

	#this is a simple wrapper for an R script so easier to run from amptk menu
	tmpdir = 'lulu_'+str(os.getpid())
	if not os.path.isdir(tmpdir):
		os.makedirs(tmpdir)

	#generate the match list using the minimum match pident
	match_file = os.path.join(tmpdir, 'match_file.txt')
	amptklib.log.info("Loading {:,} OTUs".format(amptklib.countfasta(args.fasta)))
	amptklib.log.info("Generating pairwise percent identity between OTUs using VSEARCH at {:}% identity".format(args.min_match))
	cmd = ['vsearch', '--usearch_global', os.path.abspath(args.fasta), '--db', os.path.abspath(args.fasta), '--self', '--id', str(args.min_match / 100), '--iddef', '1', '--userout', match_file, '--userfields', 'query+target+id', '--maxaccepts', '0', '--query_cov', '.9', '--maxhits', '10']
	amptklib.runSubprocess(cmd, amptklib.log)

	#now run LULU in R
	LULU_log = os.path.join(tmpdir, 'LULU-R.log')
	lulu_otu_table = base + '.lulu.otu_table.txt'
	dropList = os.path.join(tmpdir, 'droplist.txt')
	MapData = base + '.lulu.otu-map.txt'
	amptklib.log.info("Running LULU algorithm")
	cmd = ['Rscript', '--vanilla', luluScript, os.path.abspath(args.otu_table), os.path.abspath(match_file), args.min_ratio_type, str(args.min_ratio), str(args.min_match), str(args.min_relative_cooccurence / 100), lulu_otu_table, dropList, MapData]
	amptklib.runSubprocess4(cmd, amptklib.log, LULU_log)

	#get updated OTUs
	remove = []
	with open(dropList, 'rU') as dropped:
		for line in dropped:
			remove.append(line.rstrip())
	lulu_otus = base + '.lulu.otus.fa'
	with open(lulu_otus, 'w') as output:
		with open(args.fasta, 'rU') as infasta:
			for record in SeqIO.parse(infasta, 'fasta'):
				if not record.id in remove:
					output.write('>%s\n%s\n' % (record.id, record.seq))
	amptklib.log.info("LULU has merged {:,} OTUs, output data contains {:,} OTUs".format(len(remove), amptklib.countfasta(lulu_otus)))
	amptklib.log.info("LULU OTU table post processing finished\n\
----------------------------------\n\
OTU table:  {:}\n\
OTU FASTA:  {:}\n\
LULU map:   {:}\n\
----------------------------------".format(lulu_otu_table,lulu_otus, MapData))
	if 'win32' in sys.platform:
		print("\nExample of next cmd: amptk taxonomy -f %s -i %s -m mapping_file.txt -d ITS2\n" % (lulu_otus, lulu_otu_table))
	else:
		print(colr.WARN + "\nExample of next cmd:" + colr.END + " amptk taxonomy -f %s -i %s -m mapping_file.txt -d ITS2\n" % (lulu_otus, lulu_otu_table))
	if not args.debug:
		if os.path.isdir(tmpdir):
			shutil.rmtree(tmpdir)
	print("-------------------------------------------------------")

if __name__ == "__main__":
	main(args)