#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import tarfile
import gzip
import shutil
from amptk import amptklib
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

URL = { 'ITS': 'https://osf.io/pbtyh/download?version=6',
        '16S': 'https://osf.io/m7v5q/download?version=3', 
        'LSU': 'https://osf.io/sqn5r/download?version=3', 
        'COI': 'https://osf.io/pax79/download?version=4' }

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50) 
	
def main(args):
	parser=argparse.ArgumentParser(prog='amptk-install.py',
		description='''Script to download preformatted databases''',
		epilog="""Written by Jon Palmer (2019) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)
	parser.add_argument('-i','--input', nargs='+', required=True, choices=['ITS', '16S', 'LSU', 'COI'], help='Install Databases')
	parser.add_argument('-f','--force', action='store_true', help='Overwrite existing databases')
	args=parser.parse_args(args)
	
	parentdir = os.path.join(os.path.dirname(amptklib.__file__))
	
	for x in args.input:
		udbfile = os.path.join(parentdir, 'DB', x+'.udb')
		if os.path.isfile(udbfile):
			if not args.force:
				print("A formated database was found, to overwrite use '--force'. You can add more custom databases by using the `amptk database` command.")
				sys.exit(1)
		#download
		if not x in URL:
			if args.force:
				continue
			print("%s not valid, choices are ITS, 16S, LSU, COI" % x)
			sys.exit(1)
		print("Downloading %s pre-formatted database" % x)
		address = URL.get(x)
		amptklib.download(address, x+'.amptk.tar.gz')
		tfile = tarfile.open(x+'.amptk.tar.gz', 'r:gz')
		tfile.extractall(x)
		for file in os.listdir(x):
			shutil.move(os.path.join(x,file), os.path.join(parentdir, 'DB', file))
		shutil.rmtree(x)
		os.remove(x+'.amptk.tar.gz')
		print('Extracting FASTA files for {:}'.format(x))
		extracted = os.path.join(parentdir, 'DB', x+'.extracted.fa')
		cmd = ['vsearch', '--udb2fasta', udbfile, '--output', extracted]
		amptklib.runSubprocess5(cmd)
		print("{:} taxonomy database installed to {:}".format(x, os.path.join(parentdir, 'DB')))

if __name__ == "__main__":
	main()