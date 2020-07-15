#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import shutil
import glob
import argparse
import itertools
import edlib
import pyfastx
import multiprocessing
from natsort import natsorted
import amptk.amptklib as lib


class col(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'


def findFwdPrimer(primer, sequence, mismatch, equalities):
    # search for match
    align = edlib.align(primer, sequence, mode="HW",
                        k=int(mismatch), additionalEqualities=equalities)
    if align["editDistance"] >= 0:  # we found a hit
        TrimPos = align["locations"][0][1]+1
        return(TrimPos, align['editDistance'])
    else:  # return position will be False if not found
        return (False, False)


def findRevPrimer(primer, sequence, mismatch, equalities):
    # search for match
    align = edlib.align(primer, sequence, mode="HW",
                        task="locations", k=int(mismatch),
                        additionalEqualities=equalities)
    if align["editDistance"] >= 0:  # we found a hit
        TrimPos = align["locations"][0][0]
        return(TrimPos, align['editDistance'])
    else:  # return position will be False if not found
        return (False, False)


def calcOrientation(tup, length):
    # input is a tuple of tuples (loc/False,mismatch/False)
    def _count(inputTup):
        num = 0
        for x in inputTup:
            if x[0]:
                num += 1
        return num

    def _sumValues(inputTup):
        mismatch = 0
        for x in inputTup:
            if x[1]:
                mismatch += x[1]
        return mismatch
    forward = _count([tup[0], tup[2], tup[4]])
    reverse = _count([tup[1], tup[3], tup[5]])
    forCutTup = (False, False)
    revCutTup = (False, False)
    if forward == 0 and reverse == 0:  # means no primer was found
        return False, False
    elif forward == reverse:  # now need to look at mismatches
        forMatch = _sumValues([tup[0], tup[2], tup[4]])
        revMatch = _sumValues([tup[1], tup[3], tup[5]])
        if forMatch < revMatch:
            return True, [tup[2][0], tup[4][0]]
        elif revMatch < forMatch:
            if tup[3][0]:
                forCutTup = (length-tup[3][0], tup[3][1])
            if tup[5][0]:
                revCutTup = (length-tup[5][0], tup[5][1])
            return False, [forCutTup[0], revCutTup[0]]
        else:  # cant determine proper orientation, drop
            return False, False
    if forward > reverse:
        return True, [tup[2][0], tup[4][0]]
    elif reverse > forward:
        if tup[3][0]:
            forCutTup = (length-tup[3][0], tup[3][1])
        if tup[5][0]:
            revCutTup = (length-tup[5][0], tup[5][1])
        return False, [forCutTup[0], revCutTup[0]]


def fasta2dict(input):
    lookup = {}
    for i, (name, seq) in enumerate(pyfastx.Fasta(input, build_index=False)):
        lookup[i] = name
    return lookup


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: processBAM(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))


def processBAM(bam, args=False):
    def _reverse_string(s):
        return s[::-1]

    basename = os.path.basename(bam).split('.bam')[0]
    fastqOUT = os.path.join(args.out, basename+'.oriented.fq')
    StatsOut = os.path.join(args.out, basename+'.stats')
    # take care of some lookup data
    bcLookup = fasta2dict(args.barcodes)
    if args.int_primer:
        if args.int_primer == 'ITS4':
            intPrimer = lib.RevComp(lib.primer_db[args.int_primer])
        else:
            if args.int_primer in lib.primer_db:
                intPrimer = lib.primer_db[args.int_primer]
            else:
                intPrimer = args.int_primer
    else:
        intPrimer = False
    Total = 0
    NoPrimer = 0
    TooShort = 0
    BadParse = 0
    ValidSeqs = 0
    with open(StatsOut, 'w') as stats:
        with open(fastqOUT, 'w') as outfile:
            for aln in lib.execute(['samtools', 'view', bam]):
                Total += 1
                cols = aln.split('\t')
                seq = cols[9]
                qual = cols[10]
                tags = cols[11:]
                zm, rq, bc, bq = (None,)*4
                #RG:Z:59204c41  rq:f:0.99996	zm:i:66165	qs:i:17	qe:i:1361	bc:B:S,0,0	bq:i:96
                for x in tags:
                    if x.startswith('zm:'):
                        zm = x.split(':')[-1]
                    elif x.startswith('bc:'):
                        bc = x.split(':')[-1]
                        if ',' in bc:
                            code, forbc, revbc = bc.split(',')
                            if forbc == revbc:
                                bc = int(forbc)
                                bcName = bcLookup.get(bc)
                            else:
                                bc = (forbc, revbc)
                                bcName = '{}-{}'.format(bcLookup.get(forbc),
                                                        bcLookup.get(revbc))
                        else:
                            badParse += 1
                            continue
                    elif x.startswith('rq:'):
                        rq = x.split(':')[-1]
                    elif x.startswith('bq:'):
                        bq = x.split(':')[-1]
                if any(elem is None for elem in [zm, rq, bc, bq]):
                    BadParse += 1
                    continue
                newHeader = '{};barcodelabel={};bq={};rq={}'.format(zm, bcName,
                                                                     bq, rq)
                # now lets look for Primers
                # lets be safe here, not worry about speed and look for all
                if intPrimer:
                    ifPos = findRevPrimer(intPrimer, seq,
                                          args.mismatch, lib.degenNuc)
                    irPos = findRevPrimer(lib.RevComp(intPrimer), seq,
                                          args.mismatch, lib.degenNuc)
                else:
                    ifPos = (False, False)
                    irPos = (False, False)
                ffPos = findFwdPrimer(args.fwd_primer, seq, args.mismatch,
                                      lib.degenNuc)
                frPos = findRevPrimer(lib.RevComp(args.fwd_primer),
                                      seq, args.mismatch, lib.degenNuc)
                rfPos = findRevPrimer(lib.RevComp(args.rev_primer),
                                      seq, args.mismatch, lib.degenNuc)
                rrPos = findFwdPrimer(args.rev_primer, seq, args.mismatch,
                                      lib.degenNuc)
                orientation, primerMatch = calcOrientation((ifPos, irPos,
                                                            ffPos, frPos,
                                                            rfPos, rrPos),
                                                           len(seq))
                if not primerMatch:
                    NoPrimer += 1
                    continue
                if not primerMatch[0]:
                    forPos = 0
                else:
                    forPos = primerMatch[0]
                if not primerMatch[1]:
                    revPos = len(seq)
                else:
                    revPos = primerMatch[1]
                if (revPos - forPos) < args.min_len:
                    TooShort += 1
                    continue
                # now write output of any remaining
                ValidSeqs += 1
                newHeader = '{:};primers=[{:},{:}]'.format(newHeader,
                                                           primerMatch[0],
                                                           primerMatch[1])
                if orientation:
                    newHeader = '{:};orient={:}'.format(newHeader, 'plus')
                    outfile.write('@{:}\n{:}\n+\n{:}\n'.format(
                        newHeader.strip(), seq[forPos:revPos],
                        qual[forPos:revPos]))
                else:
                    newHeader = '{:};orient={:}'.format(newHeader, 'minus')
                    outfile.write('@{:}\n{:}\n+\n{:}\n'.format(
                        newHeader.strip(),
                        lib.RevComp(seq)[forPos:revPos],
                        _reverse_string(qual)[forPos:revPos]))
        stats.write('{},{},{},{},{}\n'.format(Total, BadParse, NoPrimer,
                                              TooShort, ValidSeqs))


def main(args):
    parser = argparse.ArgumentParser(
        prog='amptk-process_pacbio.py',
        description='''Script to process pacbio CCS amplicon data''',
        epilog="""Written by Jon Palmer (2020) nextgenusfs@gmail.com""")
    parser.add_argument('-i', '--bam', required=True,
                        help='Input directory containing Lima BAM files')
    parser.add_argument('-o', '--out', required=True, help='Output basename')
    parser.add_argument('-b', '--barcodes', required=True,
                        help='Barcodes FASTA file')
    parser.add_argument('-f', '--fwd_primer', default='CTTGGTCATTTAGAGGAAGTAA',
                        help='Forward Primer')
    parser.add_argument('-r', '--rev_primer', default='CCGTGTTTCAAGACGGG',
                        help='Reverse Primer')
    parser.add_argument('--int_primer', help='Internal primer')
    parser.add_argument('--primer_mismatch', dest='mismatch', default=3,
                        type=int, help='Number of mis-matches in primer')
    parser.add_argument('-m', '--min_len', default=800, type=int,
                        help='Minimum read length to keep')
    parser.add_argument('--cpus', type=int,
                        help="Number of CPUs. Default: auto")
    args = parser.parse_args(args)

    if not args.cpus:
        args.cpus = multiprocessing.cpu_count()

    # remove logfile if exists
    log_name = '{}.amptk-pacbio.log'.format(args.out)
    if os.path.isfile(log_name):
        os.remove(log_name)

    lib.setupLogging(log_name)
    cmd_args = " ".join(sys.argv)+'\n'
    lib.log.debug(cmd_args)
    print("-------------------------------------------------------")
    lib.SystemInfo()
    lib.versionDependencyChecks('usearch9')

    #create directory and check for existing logfile
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # parse the input, if directories then traverse
    lib.log.info('Expecting PacBio CCS demultiplexed BAM files (ie from ccs-->lima) in {}'.format(args.bam))
    data = []
    for f in os.listdir(args.bam):
        if f.endswith('.bam'):
            data.append(os.path.join(args.bam, f))
    lib.log.info('Found {:,} BAM files in {}'.format(len(data), args.bam))
    if args.fwd_primer in lib.primer_db:
        args.fwd_primer = lib.primer_db[args.fwd_primer]
    if args.rev_primer in lib.primer_db:
        args.rev_primer = lib.primer_db[args.rev_primer]
    if args.int_primer:
        if args.int_primer == 'ITS4':
            intPrimer = lib.RevComp(lib.primer_db[args.int_primer])
        else:
            if args.int_primer in lib.primer_db:
                intPrimer = lib.primer_db[args.int_primer]
            else:
                intPrimer = args.int_primer
        lib.log.info(
            'Reorienting reads and trimming primers; For Primer: {} Rev Primer: {} Int Primer: {}'.format(
                args.fwd_primer, args.rev_primer, intPrimer))
    else:
        lib.log.info(
            'Reorienting reads and trimming primers; For Primer: {} Rev Primer: {}'.format(
                args.fwd_primer, args.rev_primer))

    lib.runMultiProgress(safe_run, data, args.cpus, args=args)

    #parse the stats
    finalstats = [0,0,0,0,0]
    for file in os.listdir(args.out):
        if file.endswith('.stats'):
            with open(os.path.join(args.out, file), 'r') as statsfile:
                line = statsfile.readline()
                line = line.rstrip()
                newstats = line.split(',')
                newstats = [int(i) for i in newstats]
                for x, num in enumerate(newstats):
                    finalstats[x] += num
    lib.log.info('{:,} total reads'.format(finalstats[0]))
    if finalstats[2] > 0:
        lib.log.info('{:,} reads unable to parse properly'.format(finalstats[1]))
    lib.log.info('{:,} reads with no primer match'.format(finalstats[2]))
    lib.log.info('{:,} reads too short (<{} bp)'.format(finalstats[3], args.min_len))
    lib.log.info('{:,} valid output reads'.format(finalstats[4]))

    # now combine the results of properfied data
    catDemux = args.out + '.demux.fq'
    with open(catDemux, 'w') as outfile:
        for filename in glob.glob(os.path.join(args.out, '*.oriented.fq')):
            if filename == catDemux:
                continue
            with open(filename, 'r') as readfile:
                shutil.copyfileobj(readfile, outfile)

    #now loop through data and find barcoded samples, counting each.....
    BarcodeCount = {}
    with open(catDemux, 'r') as input:
        header = itertools.islice(input, 0, None, 4)
        for line in header:
            ID = line.split("label=", 1)[-1].split(";")[0]
            if ID not in BarcodeCount:
                BarcodeCount[ID] = 1
            else:
                BarcodeCount[ID] += 1

    #now let's count the barcodes found and count the number of times they are found.
    barcode_counts = "%22s:  %s" % ('Sample', 'Count')
    for k, v in natsorted(list(BarcodeCount.items()), key=lambda k_v: k_v[1], reverse=True):
        barcode_counts += "\n%22s:  %s" % (k, str(BarcodeCount[k]))
    lib.log.info("Found %i barcoded samples\n%s" % (len(BarcodeCount), barcode_counts))

    # gzip final and clean up
    FinalDemux = catDemux+'.gz'
    lib.Fzip(catDemux, FinalDemux, args.cpus)
    lib.removefile(catDemux)
    shutil.rmtree(args.out)
    #get file size
    filesize = os.path.getsize(FinalDemux)
    readablesize = lib.convertSize(filesize)
    lib.log.info("Output file:  %s (%s)" % (FinalDemux, readablesize))

    print("-------------------------------------------------------")
    if 'darwin' in sys.platform:
        print(col.WARN + "\nExample of next cmd: " + col.END + "amptk pb-dada2 -i %s -o out\n" % (FinalDemux))
    else:
        print("\nExample of next cmd: amptk pb-dada2 -i %s -o out\n" % (FinalDemux))

if __name__ == "__main__":
    main(sys.argv[1:])
