#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import logging
import shutil
import subprocess
from amptk import amptklib


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-stats.py',
		description='''Script takes BIOM as input and runs basic summary stats''',
		epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--biom', required=True, help='Input BIOM file (OTU table + metadata)')
	parser.add_argument('-t','--tree', required=True, help='Phylogentic tree from AMPtk taxonomy')
	parser.add_argument('-o','--out', default='amptk_stats', help='Output folder basename')
	parser.add_argument('-d','--distance', default='raupcrick', choices=['raupcrick','bray','unifrac','wunifrac','jaccard','aitchison', 'all'], help="Distance metric")
	parser.add_argument('--indicator_species', action='store_true', help='Run indicator species analysis')
	parser.add_argument('--ignore_otus', nargs="+", help='OTUs to drop from table and run stats')
	parser.add_argument('--ord_method', default='NMDS', choices=["DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA"], help='Ordination method')
	parser.add_argument('--ord_ellipse', action='store_true', help='Add ellipses on NMDS instead of centroids & error bars')
	#parser.add_argument('-t','--treatments', nargs='+', help='treatments (metadata variables) to run, Default: all')
	args=parser.parse_args(args)
	
	parentdir = os.path.join(os.path.dirname(amptklib.__file__))

	phyloseq_nmds = os.path.join(parentdir, 'phyloseq_nmds.R')
	parse_adonis = os.path.join(parentdir, 'parse_adonis.py')
	phyloseq_nmds_indicator = os.path.join(parentdir, 'phyloseq_nmds_indicator.R')
	if args.indicator_species:
		phyloseqCMD = phyloseq_nmds_indicator
	else:
		phyloseqCMD = phyloseq_nmds

	#remove logfile if exists
	log_name = args.out + '.amptk-stats.log'
	if os.path.isfile(log_name):
		amptklib.removefile(log_name)

	amptklib.setupLogging(log_name)
	FNULL = open(os.devnull, 'w')
	cmd_args = " ".join(sys.argv)+'\n'
	amptklib.log.debug(cmd_args)
	print("-------------------------------------------------------")
	#initialize script, log system info and usearch version
	amptklib.SystemInfo()

	#check dependencies
	programs = ['Rscript']
	amptklib.CheckDependencies(programs)
	Rversions = amptklib.checkRversion()
	R_pass = '3.2.1'
	phyloseq_pass = '1.19.1'
	#check dada2 first, if good move on, otherwise issue warning
	if not amptklib.gvc(Rversions[2], phyloseq_pass):
		amptklib.log.error("R v%s; Phyloseq v%s detected, need atleast v%s" % (Rversions[0], Rversions[2], phyloseq_pass))
		amptklib.log.error("See: https://joey711.github.io/phyloseq/index.html")
		sys.exit(1)
	amptklib.log.info("R v%s; Phyloseq v%s" % (Rversions[0], Rversions[2]))

	#this is a simple wrapper for an R script so easier to run from amptk menu
	if not os.path.isdir(args.out):
		os.makedirs(args.out)

	phylolog = os.path.join(args.out, 'phyloseq-R.log')
	if args.distance == 'all':
		distances = ['raupcrick','bray','unifrac','wunifrac','jaccard']
		amptklib.log.info("Running hypothesis test using %s distance metrics on all treatments, drawing %s for each." % (','.join(distances),args.ord_method))
		for dist in distances:
			cmd = ['Rscript', '--vanilla', phyloseqCMD, os.path.abspath(args.biom), os.path.abspath(args.tree), args.out, dist, args.ord_method, str(args.ord_ellipse)]
			if args.ignore_otus:
				cmd = cmd + args.ignore_otus
			amptklib.runSubprocess3(cmd, amptklib.log, args.out, phylolog)
	else:
		amptklib.log.info("Running hypothesis test using %s distance metric on all treatments, drawing NMDS for each." % args.distance)
		cmd = ['Rscript', '--vanilla', phyloseqCMD, os.path.abspath(args.biom), os.path.abspath(args.tree), args.out, args.distance, args.ord_method, str(args.ord_ellipse)]
		if args.ignore_otus:
			cmd = cmd + args.ignore_otus
		amptklib.runSubprocess3(cmd, amptklib.log, args.out, phylolog)
	
	#parse the adonis output
	#amptklib.log.info("Parsing p-values from hyopthesis tests generated in R")
	#subprocess.call([parse_adonis, args.out])
	amptklib.log.info('HTML output files were generated for each treatment: {:}'.format(args.out))
	print("-------------------------------------------------------")
	
if __name__ == "__main__":
	main(args)
