#!/usr/bin/env python

import sys, argparse, os, inspect, shutil, glob, multiprocessing
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import FastaIterator
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.amptklib as amptklib
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.edlib as edlib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='amptk-strip_primer.py',
    description='''Script removes primers from FASTQ files''',
    epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='Input FASTQ file')
parser.add_argument('--reverse', help='Input FASTQ file')
parser.add_argument('-f','--fwd_primer', required=True, help='Forward Primer')
parser.add_argument('-r','--rev_primer', help='Reverse Primer')
parser.add_argument('-o','--out', required=True, help='Output basename')
parser.add_argument('--primer_mismatch', default=2, type=int, help='Number of mis-matches in primer')
args=parser.parse_args()

def primerStrip(file):
    base = os.path.basename(file).split('.')[0]
    goodseq = os.path.join(tmpdir, base+'.good')
    badseq = os.path.join(tmpdir, base+'.bad')
    with open(goodseq, 'w') as good:
        with open(badseq, 'w') as bad:
            for title, seq, qual in FastqGeneralIterator(open(file)):
                foralign = edlib.align(args.fwd_primer, seq, mode="HW", k=args.primer_mismatch)
                if foralign["editDistance"] >= 0:
                    ForCutPos = foralign["locations"][0][1]+1
                    Seq = seq[ForCutPos:]
                    Qual = qual[ForCutPos:]
                    #align reverse
                    revalign = edlib.align(RevPrimer, Seq, mode="HW", task="locations", k=args.primer_mismatch)
                    if revalign["editDistance"] >= 0:
                        RevCutPos = revalign["locations"][0][0]
                        Seq = Seq[:RevCutPos]
                        Qual = Qual[:RevCutPos]                              
                    good.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))
                else:
                    bad.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))                   

def revprimerStrip(file):
    base = os.path.basename(file).split('.')[0]
    goodseq = os.path.join(tmpdir, base+'.good')
    badseq = os.path.join(tmpdir, base+'.bad')
    with open(goodseq, 'w') as good:
        with open(badseq, 'w') as bad:
            for title, seq, qual in FastqGeneralIterator(open(file)):
                foralign = edlib.align(args.rev_primer, seq, mode="HW", k=args.primer_mismatch)
                if foralign["editDistance"] >= 0:
                    ForCutPos = foralign["locations"][0][1]+1
                    Seq = seq[ForCutPos:]
                    Qual = qual[ForCutPos:]
                    #align reverse
                    revalign = edlib.align(RevForPrimer, Seq, mode="HW", task="locations", k=args.primer_mismatch)
                    if revalign["editDistance"] >= 0:
                        RevCutPos = revalign["locations"][0][0]
                        Seq = Seq[:RevCutPos]
                        Qual = Qual[:RevCutPos]                              
                    good.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))
                else:
                    bad.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))                 


def splitter(inputfile, tempdir):
    '''
    this function splits a fastq file into equal parts into a temporary directory
    and then returns the file names in a list
    '''
    total = amptklib.countfastq(inputfile)
    #split the input FASTQ file into chunks to process
    with open(inputfile, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fastq')
        chunks = total / (4*cpus)+1
        #divide into chunks, store in tmp file
        for i, batch in enumerate(amptklib.batch_iterator(SeqRecords, chunks)) :
            filename = "chunk_%i.fq" % (i+1)
            tmpout = os.path.join(tempdir, filename)
            with open(tmpout, 'w') as handle:
                SeqIO.write(batch, handle, "fastq")
    #now get file list from tmp folder
    file_list = []
    for file in os.listdir(tempdir):
        if file.endswith(".fq"):
            file = os.path.join(tempdir, file)
            file_list.append(file)
    return file_list

def combiner(dir, ending, output):
    '''
    this function takes a directory, and finds all files that have ending (.good)
    and concantenates them all into an output file
    '''
    #then concatenation good and bad then delete tmpdir
    with open(output, 'wb') as outfile:
        for filename in glob.glob(os.path.join(dir,'*'+ending)):
            if filename == output:
                continue
            with open(filename, 'rU') as readfile:
                shutil.copyfileobj(readfile, outfile)

#get cpus
cpus = multiprocessing.cpu_count()
print("-------------------------------------------------------")
#now run primer strip
if args.reverse:
    if not args.rev_primer:
        print("ERROR: if reverse reads passed you must provide -r,--rev_primer")
        sys.exit(1)
    #reverse comp primers for search
    RevPrimer = revcomp_lib.RevComp(args.rev_primer)
    RevForPrimer = revcomp_lib.RevComp(args.fwd_primer)
    #setup tmpdir
    tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    #first run forwards
    print("Splitting forward reads into buckets")
    GoodFor = args.out + '.fwd.stripped.fq'
    BadFor = args.out + '.fwd.no_primer.fq'
    filelist = splitter(args.input, tmpdir)
    print("Stripping primers from forward reads")
    amptklib.runMultiProgress(primerStrip, filelist, cpus)
    combiner(tmpdir, ".good", GoodFor)
    combiner(tmpdir, ".bad", BadFor)
    shutil.rmtree(tmpdir)

    #now run reverse
    print("Splitting reverse reads into buckets")
    tmpdir = args.out.split('.')[0]+'_'+str(os.getpid()+1)
    os.makedirs(tmpdir)
    GoodRev = args.out + '.rev.stripped.fq'
    BadRev = args.out + '.rev.no_primer.fq'
    revfilelist = splitter(args.reverse, tmpdir)
    print("Stripping primers from reverse reads")
    amptklib.runMultiProgress(revprimerStrip, revfilelist, cpus)
    combiner(tmpdir, ".good", GoodRev)
    combiner(tmpdir, ".bad", BadRev)
    shutil.rmtree(tmpdir)

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
    print("-------------------------------------------------------")
    print("%i total reads" % total)
    print("%i primer found properly paired: %s, %s" % (passed, PEfor, PErev))
    print("%i primer found singletons: %s, %s" % (nopaired, SEfor, SErev))
    print("%i primer not found in either forward or reverse reads" % (failed))    
    amptklib.removefile(GoodFor)
    amptklib.removefile(GoodRev)
else:
    #setup tmpdir
    tmpdir = args.out.split('.')[0]+'_'+str(os.getpid())
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
     
    GoodOut = args.out +'.stripped.fq'
    BadOut = args.out +'.no_primer_found.fq'
    filelist = splitter(args.input, tmpdir)
    amptklib.runMultiProgress(primerStrip, filelist, cpus)
    combiner(tmpdir, ".good", GoodOut)
    combiner(tmpdir, ".bad", BadOut)
    shutil.rmtree(tmpdir)
    total = amptklib.countfastq(args.input)
    passed = amptklib.countfastq(GoodOut)
    failed = amptklib.countfastq(BadOut)
    print("-------------------------------------------------------")
    print("\n%i total reads" % total)
    print("%i primer found/stripped: %s" % (passed, GoodOut))
    print("%i primer not found: %s" % (failed, BadOut))
print("-------------------------------------------------------")