#!/usr/bin/env python

#This script filters results from ficus-OTU_cluster.py
#written by Jon Palmer palmer.jona at gmail dot com

import os
import argparse
import subprocess
import inspect

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-mock_filter.py',
    description='''Script inspects output of ficus-OTU_cluster.py and 
    determines useful threshold for OTU output based on a spike-in 
    mock community.''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('base', help='BaseName of ficus-OTU_cluster.py')
parser.add_argument('-m','--mock', default='ufits_mock3.fa', help='Multi-FASTA mock community')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

#open log file for usearch8 stderr redirect
log_name = 'ufits-mock.log'
if os.path.isfile(log_name):
    os.remove(log_name)
log_file = open(log_name, 'ab')

usearch = args.usearch
try:
    subprocess.call([usearch, '--version'], stdout = log_file, stderr = log_file)
except OSError:
    print "%s not found in your PATH, exiting." % usearch 
    os._exit(1)
    
