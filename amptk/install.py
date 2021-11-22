#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import tarfile
import gzip
import json
import requests
import shutil
from amptk import amptklib
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

def main(args):
    parser=argparse.ArgumentParser(prog='amptk-install.py',
        description='''Script to download preformatted databases''',
        epilog="""Written by Jon Palmer (2019) nextgenusfs@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('-i','--input', nargs='+', required=True, choices=['ITS', '16S', 'LSU', 'COI', 'PR2'], help='Install Databases')
    parser.add_argument('-f','--force', action='store_true', help='Overwrite existing databases')
    parser.add_argument('-l','--local', action='store_true', help='Use local downloads.json for links')
    args=parser.parse_args(args)

    parentdir = os.path.join(os.path.dirname(amptklib.__file__))

    # downd from github to get most recent databases
    if not args.local:
        try:
            print('Retrieving download links from GitHub Repo')
            URL = json.loads(requests.get("https://raw.githubusercontent.com/nextgenusfs/amptk/master/amptk/downloads.json").text)
        except:
            print('Unable to download links from GitHub, using funannotate version specific links')
            with open(os.path.join(os.path.dirname(__file__), 'downloads.json')) as infile:
                URL = json.load(infile)
    else:
        with open(os.path.join(os.path.dirname(__file__), 'downloads.json')) as infile:
            URL = json.load(infile)

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
            print("%s not valid, choices are ITS, 16S, LSU, COI, PR2" % x)
            sys.exit(1)
        print("Downloading %s pre-formatted database" % x)
        address = URL.get(x)
        if not os.path.isfile(x+'.amptk.tar.gz'):
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
