#!/usr/bin/env python

import sys, re, argparse, os, subprocess, shutil
from natsort import natsorted
from Bio.SeqIO.FastaIO import FastaIterator

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='bold2ufits.py',
    description='''Parse bold2utax.py output, search for primer, truncate, and construct DB''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='FASTA file from bold2utax.py')
parser.add_argument('-o','--out', required=True, help='Basename for output files')
parser.add_argument('-p','--primer', default='GGTCAACAAATCATAAAGATATTGG', help='Forward Primer Sequence')
parser.add_argument('--primer_mismatch', default='4', help='Mismatches allowed in primer')
parser.add_argument('--utax_size', default='105000', help='Number of seqs for UTAX training')
parser.add_argument('-l','--trunclen', type=int, default=200, help='Length to truncate sequences to')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

FNULL = open(os.devnull, 'w')
pid = os.getpid()
primer = 'bold2ufits_'+str(pid)+'.primer.tmp'
primermatch = 'bold2ufits_'+str(pid)+'.match.tmp'

#first create tmp primer fasta file
with open(primer, 'w') as primerout:
    primerout.write('>COIprimer\n%s\n' % args.primer)

#now search for primers using usearch
subprocess.call([args.usearch, '-search_oligodb', args.input, '-db', primer, '-strand', 'plus', '-userout', primermatch, '-userfields', 'query+target+qstrand+diffs+tlo+thi+trowdots', '-notrunclabels', '-maxdiffs', args.primer_mismatch], stdout = FNULL, stderr = FNULL)

#load results into dictionary
PrimerFound = {}
with open(primermatch, 'rU') as input:
    for line in input:
        cols = line.split('\t')
        ID = cols[0]
        ID = ID.split(';')[0]
        trunc = cols[5]
        if ID not in PrimerFound:
            PrimerFound[ID] = int(trunc)

#loop through seqs, remove primer if found, and truncate to length
truncated = 'bold2ufits_'+str(pid)+'.truncate.tmp'
with open(truncated, 'w') as output:
    for record in FastaIterator(open(args.input)):
        RecID = record.id.split(';')[0]
        if RecID in PrimerFound:
            pos = PrimerFound.get(RecID)
            Seq = record.seq[pos:args.trunclen]
        else:
            Seq = record.seq[:args.trunclen]
        output.write('>%s\n%s\n' % (record.description, Seq))

#look for exact replicates via dereplication
derep_out = 'bold2ufits_'+str(pid)+'.derep.tmp'
subprocess.call([args.usearch, '-derep_fulllength', truncated, '-uc', derep_out, '-notrunclabels'], stdout = FNULL, stderr = FNULL)
DerepRemove = []
with open(derep_out, 'rU') as dereplication:
    for line in dereplication:
        line = line.replace('\n', '')
        if line.endswith('*'):
            continue
        cols = line.split('\t')
        ID1tax = cols[8].split(';')[1]
        ID2tax = cols[9].split(';')[1]
        if ID1tax == ID2tax:
            DerepRemove.append(cols[8].split(';')[0])
            continue
        if ID1tax.count(',') > ID2tax.count(','):
            DerepRemove.append(cols[9].split(';')[0])
        else:
            DerepRemove.append(cols[8].split(';')[0])

#finally output the dereplicated/truncated sequences to output
usearch = args.out+'.all4usearch.fa'
utax = args.out+'.genus4utax.fa'
total_count = 0
with open(usearch, 'w') as output:
    with open(utax, 'w') as output2:
        for record in FastaIterator(open(truncated)):
            RecID = record.id.split(';')[0]
            if not RecID in DerepRemove: #get rid of those in derepremove
                output.write('>%s\n%s\n' % (record.description, record.seq))
                total_count += 1

#randomly subsample this dataset to 105,000 seqs, the most you can train
subprocess.call([args.usearch, '-fastx_subsample', usearch, '-sample_size', args.utax_size, '-fastaout', utax], stdout = FNULL, stderr = FNULL)
                    
print "%i total records written to %s" % (total_count, usearch)
print "%s records for UTAX training written to %s" % (args.utax_size, utax)
  
#cleanup
os.remove(primer)
os.remove(primermatch)
os.remove(derep_out)
os.remove(truncated)
