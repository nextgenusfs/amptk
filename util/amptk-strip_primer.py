#!/usr/bin/env python

import sys, argparse, os, inspect
from Bio.SeqIO.QualityIO import FastqGeneralIterator
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.amptklib as amptklib
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='amptk-strip_primer.py',
    description='''Script removes primers from FASTQ files''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input FASTQ file')
parser.add_argument('--reverse', required=True, help='Input FASTQ file')
parser.add_argument('-f','--fwd_primer', required=True, help='Forward Primer')
parser.add_argument('-r','--rev_primer', help='Reverse Primer')
parser.add_argument('-o','--out', required=True, help='Output basename')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
args=parser.parse_args()

def primerStrip(file, GoodOut, BadOut, fwdprimer, revprimer):
    PL = len(fwdprimer)
    with open(GoodOut, 'w') as good:
        with open(BadOut, 'w') as bad:
            for title, seq, qual in FastqGeneralIterator(open(file)):
                Diffs = primer.MatchPrefix(seq, fwdprimer)
                if Diffs <= args.primer_mismatch:
                    Seq = seq[PL:]
                    Qual = qual[PL:]
                    if revprimer:#now need to look for reverse primer
                        BestPosRev, BestDiffsRev = primer.BestMatch2(Seq, revcomp_lib.RevComp(revprimer), args.primer_mismatch)
                        if BestPosRev > 0:  #reverse primer was found
                            Seq = Seq[:BestPosRev]
                            Qual = Qual[:BestPosRev]                                           
                    good.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))
                else:
                    bad.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))                   

def primer2Strip(file, GoodOut, BadOut, fwdprimer, revprimer):


#now run primer strip
if args.reverse:
    if not args.rev_primer:
        print("ERROR: if reverse reads passed you must provide -r,--rev_primer")
        sys.exit(1)
    #first run forwards
    GoodFor = args.out + '.fwd.stripped.fq'
    BadFor = args.out + '.fwd.no_primer.fq'
    primerStrip(args.input, GoodFor, BadFor, args.fwd_primer, args.rev_primer)
    #now run reverse
    GoodRev = args.out + '.rev.stripped.fq'
    BadRev = args.out + '.rev.no_primer.fq'
    primerStrip(args.reverse, GoodRev, BadRev, args.rev_primer, args.fwd_primer)
    #now get bad reads into list
    singleRev = []
    singleFor = []
    for title, seq, qual in FastqGeneralIterator(open(BadFor)):
        singleRev.append(title.split(' ')[0])
    for title, seq, qual in FastqGeneralIterator(open(BadRev)):
        singleFor.append(title.split(' ')[0])
    bothfail = sorted(set(singleRev) & set(singleFor), key = singleFor.index)

    #now get PE and singletons
    PEfor = args.out +'.pe_R1.fastq'
    PErev = args.out +'.pe_R2.fastq'
    SEfor = args.out +'.fwd.singletons.fastq'
    SErev = args.out +'.rev.singletons.fastq'
    with open(PEfor, 'w') as peF:
        with open(SEfor, 'w') as seF:
            for title, seq, qual in FastqGeneralIterator(open(GoodFor)):
                ID = title.split(' ')[0]
                if not ID in singleFor:
                    peF.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                else:
                    seF.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    with open(PErev, 'w') as peR:
        with open(SErev, 'w') as seR:
            for title, seq, qual in FastqGeneralIterator(open(GoodRev)):
                ID = title.split(' ')[0]
                if not ID in singleRev:
                    peR.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                else:
                    seR.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    #do some counts of output and cleanup
    total = amptklib.countfastq(args.input)
    passed = amptklib.countfastq(PEfor)
    passedrev = amptklib.countfastq(PErev)
    if passed != passedrev:
        print("Error: forward reads %i != reverse reads %i" % (passed, passedrev))
    nopaired = len(singleFor) + len(singleRev)
    failed = len(bothfail)
    print("%i total reads" % total)
    print("%i primer found properly paired: %s, %s" % (passed, PEfor, PErev))
    print("%i primer found singletons: %s, %s" % (nopaired, SEfor, SErev))
    print("%i primer not found in either forward or reverse reads" % (failed))    
    amptklib.removefile(GoodFor)
    amptklib.removefile(GoodRev)
else:        
    GoodOut = args.out +'.stripped.fq'
    BadOut = args.out +'.no_primer_found.fq'
    primerStrip(args.input, GoodOut, BadOut, args.fwd_primer, False)
    total = amptklib.countfastq(args.input)
    passed = amptklib.countfastq(GoodOut)
    failed = amptklib.countfastq(BadOut)
    print("%i total reads" % total)
    print("%i primer found/stripped: %s" % (passed, GoodOut))
    print("%i primer not found: %s" % (failed, BadOut))
