#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import shutil
import glob
import argparse
import subprocess
import edlib
import pyfastx
import itertools
import amptk.amptklib as lib


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


def proper(input, output, minlen=600, mismatch=4, primer_required='off'):
    proper = 0
    forward = 0
    reverse = 0
    total = 0
    with open(output, 'w') as outfile:
        for header, seq, qual in pyfastx.Fastq(input, build_index=False):
            total += 1
            # parse header to get cut locations
            fPos, rPos = (False,)*2
            for y in header.split(';'):
                if y.startswith('primers='):
                    fPos, rPos = eval(y.replace('primers=', ''))
            if fPos and rPos and rPos-fPos >= minlen:
                outfile.write('@{:};\n{:}\n+\n{:}\n'.format(
                    header, seq[fPos:rPos], qual[fPos:rPos]))
                proper += 1
            elif fPos and not rPos:
                forward += 1
                if primer_required == 'off':
                    outfile.write('@{:};\n{:}\n+\n{:}\n'.format(
                        header, seq[fPos:rPos], qual[fPos:rPos]))
            elif rPos and not fPos:
                reverse += 1
                if primer_required == 'off':
                    outfile.write('@{:};\n{:}\n+\n{:}\n'.format(
                        header, seq[fPos:rPos], qual[fPos:rPos]))
            else:
                if primer_required == 'off':
                    outfile.write('@{:};\n{:}\n+\n{:}\n'.format(
                        header, seq[fPos:rPos], qual[fPos:rPos]))
    return [total, proper, forward, reverse]


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


def orient(input, output, fwdprimer='TCCGTAGGTGAACCTGCGG',
           revprimer='CCGTGTTTCAAGACGGG', internalprimer=False,
           mismatch=4, minlen=800, minqual=7, barcode=False,
           primer_required='on'):
    found = 0
    reverse = 0
    nothing = 0
    tooshort = 0
    lowquality = 0
    total = 0

    def _reverse_string(s):
        return s[::-1]

    def _ave_qual(quals):
        from math import log
        if quals:
            return -10 * log(sum([10**(q / -10) for q in quals]) / len(quals), 10)
        else:
            return None

    def _clean_header(title, quality, length):
        #@37b32633-a661-4703-9b0e-96c0d456e693 runid=7079c0cfa5bc9d05431a04a143727960694e8694 read=1044 ch=152 start_time=2020-02-09T19:32:37Z flow_cell_id=FAM95706 protocol_group_id=Feb09_Its sample_id=Feb09_Its barcode=barcode01
        info = title.split()
        barcodeID = barcode
        flowcell = False
        yacrd = ''
        for x in info:
            if x.startswith('barcode='):
                if not barcodeID:
                    barcodeID = x.replace('barcode=', '')
            elif x.startswith('flow_cell_id='):
                flowcell = x
            elif x.startswith('yacrd='):
                yacrd = x
        if not flowcell:
            flowcell = 'flow_cell_id=unknown'
        if yacrd:
            newHeader = '{:};barcodelabel={:};{:};len={:};avgqual={:.2f};{:}'.format(
                info[0], barcodeID, flowcell, length, quality, yacrd)
        else:
            newHeader = '{:};barcodelabel={:};{:};len={:};avgqual={:.2f}'.format(
                info[0], barcodeID, flowcell, length, quality)
        return newHeader

    with open(output, 'a') as outfile:
        for header, seq, qual in pyfastx.Fastq(input, build_index=False):
            total += 1
            if len(seq) < minlen:
                tooshort += 1
                continue
            phred = [ord(x)-33 for x in qual]
            avgQual = _ave_qual(phred)
            if avgQual < minqual:
                lowquality += 1
                continue
            newHeader = _clean_header(header, avgQual, len(seq))

            # lets be safe here, not worry about speed and look for all
            if internalprimer:
                ifPos = findRevPrimer(internalprimer, seq,
                                      mismatch, lib.degenNuc)
                irPos = findRevPrimer(lib.RevComp(internalprimer), seq,
                                      mismatch, lib.degenNuc)
            else:
                ifPos = (False, False)
                irPos = (False, False)
            ffPos = findFwdPrimer(fwdprimer, seq, mismatch, lib.degenNuc)
            frPos = findRevPrimer(lib.RevComp(fwdprimer),
                                  seq, mismatch, lib.degenNuc)
            rfPos = findRevPrimer(lib.RevComp(revprimer),
                                  seq, mismatch, lib.degenNuc)
            rrPos = findFwdPrimer(revprimer, seq, mismatch, lib.degenNuc)
            orientation, primerMatch = calcOrientation((ifPos,
                                                        irPos,
                                                        ffPos,
                                                        frPos,
                                                        rfPos,
                                                        rrPos),
                                                       len(seq))
            if not primerMatch:
                nothing += 1
                continue
            if not primerMatch[0]:
                forPos = 0
            else:
                forPos = primerMatch[0]
            if not primerMatch[1]:
                revPos = len(seq)
            else:
                revPos = primerMatch[1]
            if (revPos - forPos) < minlen:
                tooshort += 1
                continue
            # now write output of any remaining
            newHeader = '{:};primers=[{:},{:}]'.format(
                newHeader, primerMatch[0], primerMatch[1])
            if orientation:
                found += 1
                newHeader = '{:};orient={:}'.format(newHeader, 'plus')
                outfile.write('@{:}\n{:}\n+\n{:}\n'.format(
                    newHeader.strip(), seq[forPos:revPos],
                    qual[forPos:revPos]))
            else:
                reverse += 1
                newHeader = '{:};orient={:}'.format(newHeader, 'minus')
                outfile.write('@{:}\n{:}\n+\n{:}\n'.format(
                    newHeader.strip(),
                    lib.RevComp(seq)[forPos:revPos],
                    _reverse_string(qual)[forPos:revPos]))
    return [total, found, reverse, nothing, tooshort, lowquality]


def yacrdfilter(input, output, minlen=800, coverage=2, not_coverage=0.8):
    # function to run yacrd to clean up reads
    FNULL = open(os.devnull, 'w')
    overlaps = input+'.overlaps.paf'
    cmd = ['minimap2', '-x', 'ava-ont', '-g', '500', input, input]
    print('CMD: {:} > {:}'.format(' '.join(cmd), overlaps))
    with open(overlaps, 'w') as outfile:
        subprocess.call(cmd, stdout=outfile, stderr=FNULL)
    report = input+'.report.yacrd'
    cmd2 = ['yacrd', '--input', overlaps, '-c', str(coverage),
            '-n', str(not_coverage), '--output', report]
    print('CMD: {:}'.format(' '.join(cmd2)))
    subprocess.call(cmd2)
    # now parse results
    Results = {}
    with open(report, 'r') as infile:
        for line in infile:
            line = line.replace('\n', '')
            status, id, length, data = line.split('\t')
            if int(length) < minlen:
                continue
            Results[id] = {'status': status, 'trim': data}
    print(('Parsed {:} results from yacrd'.format(len(Results))))
    headernotfound = 0
    with open(output, 'w') as outfile:
        with open(input, 'r') as infile:
            for header, seq, qual in pyfastx.Fastq(infile, build_index=False):
                if ' ' in header:
                    title = header.split()[0]
                else:
                    title = header
                if title in Results:
                    if Results[title]['status'] == 'NotBad':
                        outfile.write('@{:} yacrd={:}\n{:}\n+\n{:}\n'.format(
                            header, Results[title]['trim'], seq, qual))
                else:
                    headernotfound += 1
    print(headernotfound)


def parseinfiles(input, mapping=False):
    Parsing = {}
    MapData = {}
    if mapping:
        with open(mapping, 'r') as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith('#') or line.startswith('\n'):
                    continue
                sample, folder, bc = line.split('\t')
                combinedname = '{}::{}'.format(folder, bc)
                if combinedname not in MapData:
                    MapData[combinedname] = sample
    for x in input:
        if os.path.isdir(x):
            # so this is expected, now need to walk through and find FASTQ
            for root, dirs, files in os.walk(x):
                for file in files:
                    if file.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                        subfolder = os.path.basename(root)
                        if 'unclassified' in subfolder:
                            continue
                        basename = '{}::{}'.format(x.strip(os.sep), subfolder)
                        filepath = os.path.abspath(os.path.join(root, file))
                        if basename in MapData:
                            barcodelabel = MapData[basename]
                        elif mapping:
                            continue
                        else:
                            barcodelabel = basename.replace('::', '-')
                        if basename not in Parsing:
                            Parsing[basename] = {'label': barcodelabel,
                                                 'files': [filepath]}
                        else:
                            Parsing[basename]['files'].append(filepath)
        elif os.path.isfile(x) and '.f' in x:
            basename = os.path.basename(x).rsplit('.f', 1)[0]
            Parsing[basename] = [os.path.abspath(x)]
        else:
            print('ERROR: cannot parse {}'.format(x))
    return Parsing


def main(args):
    parser = argparse.ArgumentParser(
        prog='amptk-process_ont.py',
        description='''Script to process oxford nanopore amplicon data''',
        epilog="""Written by Jon Palmer (2020) nextgenusfs@gmail.com""")
    parser.add_argument('-i', '--fastq', dest='fastq', nargs='+',
                        required=True, help='Input directory or fastq files')
    parser.add_argument('-o', '--out', required=True, help='Output basename')
    parser.add_argument('-b', '--barcode-map', dest='barcode_map',
                        help='Add barcode label')
    parser.add_argument('-f', '--fwd_primer',  default='TCCGTAGGTGAACCTGCGG',
                        help='Forward Primer')
    parser.add_argument('-r', '--rev_primer', default='CCGTGTTTCAAGACGGG',
                        help='Reverse Primer')
    parser.add_argument('--int_primer', help='Internal primer')
    parser.add_argument('--primer_mismatch', default=3, type=int,
                        help='Number of mis-matches in primer')
    parser.add_argument('-m', '--min_len', default=800, type=int,
                        help='Minimum read length to keep')
    parser.add_argument('-q', '--min_qual', default=7, type=int,
                        help='Minimum read quality to keep')
    parser.add_argument('--require_primer', default='on',
                        choices=['on', 'off'],
                        help='Require Fwd primer to be present')
    parser.add_argument('--cpus', type=int, default=1,
                        help="Number of CPUs. Default: auto")
    args = parser.parse_args(args)

    # remove logfile if exists
    log_name = '{}.amptk-ont.log'.format(args.out)
    if os.path.isfile(log_name):
        os.remove(log_name)

    lib.setupLogging(log_name)
    cmd_args = " ".join(sys.argv)+'\n'
    lib.log.debug(cmd_args)
    print("-------------------------------------------------------")
    lib.SystemInfo()

    # figure out which primers to use
    intPrimer = False
    if args.int_primer == 'ITS4':
        intPrimer = lib.RevComp(lib.primer_db[args.int_primer])
    else:
        if args.int_primer in lib.primer_db:
            intPrimer = lib.primer_db[args.int_primer]
        else:
            intPrimer = args.int_primer
    if args.fwd_primer in lib.primer_db:
        args.fwd_primer = lib.primer_db[args.fwd_primer]
    if args.rev_primer in lib.primer_db:
        args.rev_primer = lib.primer_db[args.rev_primer]

    # setup tmpdir storage
    global tmpdir
    tmpdir = 'proc_ont_{}'.format(os.getpid())
    os.makedirs(tmpdir)

    # parse the input, if directories then traverse
    if os.path.isdir(args.fastq):
        data = parseinfiles(args.fastq, mapping=args.barcode_map)
        step1Stats = []
        for k, v in data.items():
            barcode = v['label']
            run, guppybarcode = k.split('::')
            lib.log.info('Working on {:,} FASTQ files from {} {} --> {}'.format(
                len(v['files']), run, guppybarcode, barcode))
            for reads in v['files']:
                oriented = os.path.join(tmpdir, '{}.{}.oriented.fq'.format(
                    run, guppybarcode))
                lib.log.debug('Parsing: {}'.format(reads))
                orientStats = orient(reads,
                                    oriented,
                                    fwdprimer=args.fwd_primer,
                                    revprimer=args.rev_primer,
                                    internalprimer=intPrimer,
                                    mismatch=args.primer_mismatch,
                                    minlen=args.min_len,
                                    minqual=args.min_qual,
                                    barcode=barcode)
                lib.log.debug('Demuxed: {}'.format(orientStats))
                step1Stats.append(orientStats)
    elif os.path.isfile(args.fastq):  # here passed single file so need to orient and find barcodes

    combinedStats = [sum(i) for i in zip(*step1Stats)]
    lib.log.info('{:,} of {:,} reads passed'.format(
        combinedStats[1]+combinedStats[2], combinedStats[0]))
    lib.log.info('  {:,} could not orient'.format(combinedStats[3]))
    lib.log.info('  {:,} too short'.format(combinedStats[4]))
    lib.log.info('  {:,} low quality'.format(combinedStats[5]))
    # [total, found, reverse, nothing, tooshort, lowquality]
    # now combine the results of properfied data
    catDemux = args.out + '.demux.fq'
    with open(catDemux, 'w') as outfile:
        for filename in glob.glob(os.path.join(tmpdir, '*.oriented.fq')):
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

    FinalDemux = catDemux+'.gz'
    lib.Fzip(catDemux, FinalDemux, args.cpus)
    lib.removefile(catDemux)
    shutil.rmtree(tmpdir)

    #get file size
    filesize = os.path.getsize(FinalDemux)
    readablesize = lib.convertSize(filesize)
    lib.log.info("Output file:  %s (%s)" % (FinalDemux, readablesize))

    print("-------------------------------------------------------")
    if 'darwin' in sys.platform:
        print(col.WARN + "\nExample of next cmd: " + col.END + "amptk ont-extend -i %s -f seeds.fa\n" % (FinalDemux))
    else:
        print("\nExample of next cmd: amptk ont-extend -i %s -f seeds.fa\n" % (FinalDemux))


if __name__ == "__main__":
    main(sys.argv[1:])
