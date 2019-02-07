#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import range
import sys
import re
import argparse
import os
import subprocess
import shutil
import edlib
import io
from natsort import natsorted
from Bio.SeqIO.FastaIO import FastaIterator
from Bio import SeqIO
from amptk import amptklib
from future.utils import native

def dereplicate(input, output):
    seqs = {}
    with open(output, 'w') as out:
        with open(input, 'rU') as in_file:
            for rec in SeqIO.parse(in_file, 'fasta'):
                sequence = str(rec.seq)
                if not sequence in seqs:
                    seqs[sequence] = rec.description
                else:
                    #check length of taxonomy string, keep one with more tax info
                    newHeader = rec.description.split(';tax=')
                    oldHeader = seqs.get(sequence).split(';tax=')
                    newTax = newHeader[-1].split(',')
                    oldTax = oldHeader[-1].split(',')
                    newID = newHeader[0]
                    oldID = oldHeader[0]
                    newTaxLen = len(newTax)
                    oldTaxLen = len(oldTax)
                    if newTaxLen > oldTaxLen:
                        seqs[sequence] = rec.description
                    elif newTaxLen == oldTaxLen:
                        if newTaxLen == 1: #this means we have only a single level of taxonomy, so just move to next record
                            continue
                        if newTax[-1] == oldTax[-1]: #so taxonomy is the same, keep current value in dict and move to next record
                            continue
                        else: #loop backwards through tax string find last common ancestor
                            #amptklib.log.debug("ERROR: %s and %s have identical sequences, but taxonomy doesn't agree" % (','.join(oldTax), ','.join(newTax)))
                            lca = 0
                            for num in range(1,newTaxLen+1):
                                if newTax[-num] == oldTax[-num]:
                                    lca = num-1
                                    break
                            consensusTax = oldID+';tax='+','.join(oldTax[:-lca])
                            #amptklib.log.debug("setting taxonomy to %s" % (consensusTax))
                            seqs[sequence] = consensusTax
        #now write to file
        seenID = []
        seenTax = []
        for key,value in seqs.items():
            ID = value.split(';tax=')[0]
            tax = value.split(';tax=')[-1]
            taxCount = tax.count(',')
            if ID in seenID: #don't do duplicate IDs
                continue
            if taxCount > 5: #species level, not have more than 1 seq per "species"
                if ID in seenTax:
                    continue
            out.write('>'+value+'\n'+key+'\n')
            seenTax.append(tax)
            seenID.append(ID)
            
def main():
	#setup menu with argparse
	class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
		def __init__(self,prog):
			super(MyFormatter,self).__init__(prog,max_help_position=48)
	parser=argparse.ArgumentParser(prog='bold2amptk.py',
		description='''Parse bold2utax.py output, search for primer, truncate, and construct DB''',
		epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
		formatter_class = MyFormatter)
	parser.add_argument('-i','--input', required=True, help='FASTA file from bold2utax.py')
	parser.add_argument('-o','--out', required=True, help='Basename for output files')
	parser.add_argument('-f','--fwdprimer', default='GGTCAACAAATCATAAAGATATTGG', help='Forward Primer Sequence')
	parser.add_argument('-r','--revprimer', default='GGRGGRTASACSGTTCASCCSGTSCC', help='Forward Primer Sequence')
	parser.add_argument('--primer_mismatch', default=4, type=int, help='Mismatches allowed in primer')
	parser.add_argument('--utax_size', default=90000, type=int, help='Number of seqs for UTAX training')
	parser.add_argument('-m','--minlen', type=int, default=200, help='Minimum length')
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
	args=parser.parse_args()
				
	FNULL = open(os.devnull, 'w')
	pid = os.getpid()
	#reverse complement rev primer
	ForPrimer = args.fwdprimer
	RevPrimer = amptklib.RevComp(args.revprimer)

	print('Loading '+'{0:,}'.format(amptklib.countfasta(args.input)) + ' sequence records')
	print('Searching for forward primer: %s, and reverse primer: %s' % (ForPrimer, RevPrimer))
	print('Requiring reverse primer match with at least %i mismatches' % args.primer_mismatch)
	#loop through seqs, remove primer if found, and truncate to length
	truncated = 'bold2amptk_'+str(pid)+'.truncate.tmp'
	with open(truncated, 'w') as output:
		for record in FastaIterator(io.open(args.input)):
			Seq = str(record.seq)
			StripSeq = ''
			ForCutPos = amptklib.findFwdPrimer(native(ForPrimer), native(Seq), args.primer_mismatch, amptklib.degenNucSimple)
			RevCutPos = amptklib.findRevPrimer(native(RevPrimer), native(Seq), args.primer_mismatch, amptklib.degenNucSimple)
			if ForCutPos and RevCutPos:
				StripSeq = Seq[ForCutPos:RevCutPos]
			elif not ForCutPos and RevCutPos:
				StripSeq = Seq[:RevCutPos]
			if len(StripSeq) >= args.minlen:
				output.write('>%s\n%s\n' % (record.description, amptklib.softwrap(StripSeq)))


	#look for exact replicates via dereplication
	trunctotal = amptklib.countfasta(truncated)
	derep_out = 'bold2amptk_'+str(pid)+'.derep.tmp'
	usearch = args.out+'.all4usearch.fa'
	utax = args.out+'.genus4utax.fa'
	print('{0:,}'.format(trunctotal)+ ' seqs passed -> now dereplicating')
	dereplicate(truncated, usearch)
	dereptotal = amptklib.countfasta(usearch)

	#randomly subsample this dataset to make sure can train with UTAX
	if dereptotal > args.utax_size:
		subprocess.call([args.usearch, '-fastx_subsample', usearch, '-sample_size', str(args.utax_size), '-fastaout', utax], stdout = FNULL, stderr = FNULL)
	else:
		shutil.copyfile(usearch, utax) 
	utaxCount = amptklib.countfasta(utax)       
	print("{:,} unique records written to {:}".format(dereptotal, usearch))
	print("{:,} records for UTAX training written to {:}".format(utaxCount, utax))
  
	#cleanup
	os.remove(truncated)
	
if __name__ == "__main__":
	main()
