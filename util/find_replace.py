#!/usr/bin/env python

#script to do quick find and replace using dictionary

import argparse
import re

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='find_replace.py',
    description='''Script finds and replaces using dictionary lookup''',
    epilog="""Written by Jon Palmer (2015)  palmer.jona@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', dest="input", required=True, help='input File (Required)')
parser.add_argument('-d','--dictionary', dest="dict", required=True, help='Dictionary file (2 columns)')
parser.add_argument('-o','--out', default='out', help='Output file')
parser.add_argument('--verbose', action="store_true")
args=parser.parse_args()
  
def replace_all(text, wordDict):
    for key in wordDict:
        text = text.replace(key, wordDict[key])
    return text

#load dict file
dictionary = {}
in_file = open(args.dict, "rb")
for line in in_file:
    line = line.split("\t")
    line[-1] = line[-1].strip("\n")
    if line:
        dictionary[line[0]]=line[1]
    else:
        continue
if args.verbose:
    print "Dictionary"
    print dictionary

#load input file and write new one
output_file = open(args.out, "w")
with open(args.input, "r") as myinput:
    for line in myinput:
        new_line = replace_all(line, dictionary)
        output_file.write(str(new_line))

output_file.close()
