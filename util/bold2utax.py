#!/usr/bin/env python

import sys, re, argparse, os, subprocess, shutil
from natsort import natsorted
from Bio import pairwise2
from Bio.SeqIO.FastaIO import FastaIterator

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='bold2utax.py',
    description='''Parse BOLD DB TSV data dump into FASTA with UTAX compatible labels.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Bold data dump TSV format')
parser.add_argument('-o','--out', required=True, help='UTAX formated FASTA output')
parser.add_argument('--require_genbank', action='store_true', help='Require output to have GenBank Accessions')
args=parser.parse_args()

LetterToSet = {}
LetterToSet['A'] = "A"
LetterToSet['C'] = "C"
LetterToSet['G'] = "G"
LetterToSet['T'] = "T"
LetterToSet['M'] = "AC"
LetterToSet['R'] = "AG"
LetterToSet['W'] = "AT"
LetterToSet['S'] = "CG"
LetterToSet['Y'] = "CT"
LetterToSet['K'] = "GT"
LetterToSet['V'] = "ACG"
LetterToSet['H'] = "ACT"
LetterToSet['D'] = "AGT"
LetterToSet['B'] = "CGT"
LetterToSet['X'] = "GATC"
LetterToSet['N'] = "GATC"

def MatchLetter(a, b):
	global LetterToSet
	try:
		sa = LetterToSet[a.upper()]
	except:
		return False
	try:
		sb = LetterToSet[b.upper()]
	except:
		return False
	for ca in sa:
		if ca in sb:
			return True
	return False

def MatchPrefix(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	n = PrimerLength
	if L < n:
		n = L
	Diffs = 0
	for i in range(0, n):
		if not MatchLetter(Seq[i], Primer[i]):
			Diffs += 1
	return Diffs

def BestMatch(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	BestDiffs = PrimerLength
	BestPos = -1
	for Pos in range(0, L-PrimerLength+1):
		d = MatchPrefix(Seq[Pos:], Primer)
		if d < BestDiffs:
			BestDiffs = d
			BestPos = Pos
	return BestPos, BestDiffs

def TrimPrimer(Sequence, primer, mismatch):
    #find primer location
    Diffs = BestMatch(Sequence, primer)
    if Diffs[1] > int(mismatch):
        #assume that primer was trimmed off
        Seq = Sequence
    else:
        Seq = Sequence[Diffs[0]+len(primer):]
    return Seq

def align_sequences(sequence_A, sequence_B, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """
    def _calculate_identity(sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """
        sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
        matches = [sa[i] == sb[i] for i in xrange(sl)]
        seq_id = (100 * sum(matches)) / sl

        gapless_sl = sum([1 for i in xrange(sl) if (sa[i] != '-' and sb[i] != '-')])
        gap_id = (100 * sum(matches)) / gapless_sl
        return (seq_id, gap_id)

    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    alns = pairwise2.align.globalxs(sequence_A, sequence_B,
                                    gap_open, gap_extend)
    #get best alignment
    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln
    # Calculate sequence identity
    seq_id, g_seq_id = _calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id)


Total = -1
nonCOI = 0
noBIN = 0
count = 0
#create tmpdirectory
pid = os.getpid()
tmp = 'boldutax_'+str(pid)
if not os.path.exists(tmp):
    os.makedirs(tmp)

#store taxonomy in dictionary
BINtax = {}
with open(args.input, 'rU') as input:
    for line in input:
        Total += 1
        line = line.replace('\n', '')
        if line.startswith('processid'):
            header = line.split('\t')
            pid = header.index('phylum_name')
            cid = header.index('class_name')
            oid = header.index('order_name')
            fid = header.index('family_name')
            gid = header.index('genus_name')
            sid = header.index('species_name')
            seqid = header.index('nucleotides')
            boldid = header.index('sequenceID')
            gbid = header.index('genbank_accession')
            bin = header.index('bin_uri')
            idby = header.index('identification_provided_by')
            marker = header.index('marker_codes')
            continue
        #split each line at tabs
        col = line.split('\t')

        #apparently there are other genes in here, so ignore anything not COI
        if not 'COI' in col[marker]:
            nonCOI += 1
            continue

        #check for BIN, if none, then move on
        BIN = col[bin].strip()
        if BIN == '':
            noBIN += 1
            continue

        #some idiots have collector names in these places in the DB, this DB is kind of a mess, doing the best I can....
        bd = col[idby].strip()
        badnames = bd.split(' ')
        badfiltered = ['1','2','3','4','5','6','7','8','9','0']
        for y in badnames:
            if len(y) > 2:
                if y != 'Art':
                    if y != 'Eric':
                        badfiltered.append(y)
        K = 'k:Animalia'
        P = 'p:'+col[pid].strip()
        C = 'c:'+col[cid].strip()
        O = 'o:'+col[oid].strip()
        F = 'f:'+col[fid].strip()
        G = 'g:'+col[gid].strip()
        S = 's:'+col[sid].strip().replace('.', '')
        if ' sp ' in S: #remove those that have sp. in them
            S = ''
        if S.endswith(' sp'):
            S = ''
        if badfiltered:
            if any(bad in G for bad in badfiltered):
                G = ''
            if any(bad in S for bad in badfiltered):
                S = ''
        ID = col[boldid].strip()
        GB = col[gbid].strip()
        if args.require_genbank:
            if GB: #if there is a GB accession
                if 'Pending' in GB:
                    continue
                else:
                    pass
            else:
                continue
        #clean up sequence, remove any gaps, remove terminal N's
        Seq = col[seqid].replace('-', '')
        Seq = Seq.rstrip('N')
        Seq = Seq.lstrip('N')
        #if still N's in sequence, just drop it
        if 'N' in Seq:
            continue
        #get taxonomy information
        tax = []
        for i in [K,P,C,O,F,G,S]:
            if not i.endswith(':'):
                tax.append(i)
        tax_fmt = ','.join(tax)
        tax_fmt = tax_fmt.rstrip()
        if tax_fmt.endswith(','):
            tax_fmt = tax_fmt.rsplit(',',1)[0]
        if not BIN in BINtax:
            BINtax[BIN] = tax_fmt
        else:
            oldtax = BINtax.get(BIN)
            oldtax_count = oldtax.count(',')
            newtax_count = tax_fmt.count(',')
            if newtax_count > oldtax_count:
                BINtax[BIN] = tax_fmt
        #just write to fasta file first
        count += 1
        name = BIN.replace(':', '_')
        with open(os.path.join(tmp, name+'.fa'), 'ab') as output:
            if GB:
                output.write('>%s_%s;tax=%s\n%s\n' % (BIN, GB, tax_fmt, Seq))
            else:
                output.write('>%s_NA;tax=%s\n%s\n' % (BIN, tax_fmt, Seq))
        
print "%i total records processed" % Total
print "%i non COI records dropped" % nonCOI
print "%i records without a BIN dropped" % noBIN
print "%i records written to BINs" % count
print "Now looping through BINs and clustering with VSEARCH"
FNULL = open(os.devnull, 'w')
for file in os.listdir(tmp):
    base = file.replace('.fa', '')
    cluster_out = os.path.join(tmp, base+'.centroids.fa')
    subprocess.call(['vsearch', '--cluster_fast', os.path.join(tmp, file), '--id', '0.97', '--consout', cluster_out, '--notrunclabels'], stdout = FNULL, stderr = FNULL)
print "Combining consensus for each BIN"
#now grab all the centroids and combine
with open(args.out+'.tmp', 'w') as tmpout:
    for file in os.listdir(tmp):
        if file.endswith('.centroids.fa'):
            with open(os.path.join(tmp, file)) as input:
                tmpout.write(input.read()) 

print "Updating taxonomy"
#finally loop through centroids and get taxonomy from dictionary 
finalcount = 0
with open(args.out, 'w') as outputfile:
    for record in FastaIterator(open(args.out+'.tmp')):
        record.id = record.id.replace('centroid=', '')
        finalcount += 1
        fullname = record.id.split(';')[0]
        ID = fullname.split('_')[0]
        tax = BINtax.get(ID)
        outputfile.write('>%s;tax=%s\n%s\n' % (ID, tax, record.seq)) 

print "Wrote %i consensus seqs for each BIN to %s" % (finalcount, args.out)
shutil.rmtree(tmp)
os.remove(args.out+'.tmp')
