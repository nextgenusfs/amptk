#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import pyfastx
from amptk import amptklib


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter,self).__init__(prog, max_help_position=50)


def countBarcodes(file):
    # now loop through data and find barcoded samples, counting each.....
    BarcodeCount = {}
    for title, seq, qual in pyfastx.Fastq(file, build_index=False):
        if 'label=' in title:
            ID = title.split("label=", 1)[-1].split(";")[0]
        elif 'sample=' in title:
            ID = title.split("sample=", 1)[-1].split(";")[0]
        else:
            ID = title.split("=", 1)[-1].split(";")[0]
        if ID not in BarcodeCount:
            BarcodeCount[ID] = 1
        else:
            BarcodeCount[ID] += 1
    return BarcodeCount


def filter_sample(file, keep_list, output, format='fastq'):
    keep_count = 0
    total_count = 0
    with open(output, 'w') as out:
        for title, seq, qual in pyfastx.Fastq(file, build_index=False):
            total_count += 1
            if 'label=' in title:
                sample = title.split('label=', 1)[1].split(';')[0]
            elif 'sample=' in title:
                sample = title.split('sample=', 1)[1].split(';')[0]
            else:
                sample = title.split('=', 1)[1].split(';')[0]
            if sample in keep_list:
                keep_count += 1
                if format == 'fastq':
                    out.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                if format == 'fasta':
                    out.write(">%s\n%s\n" % (title, seq))
    return keep_count, total_count


def main(args):
    parser = argparse.ArgumentParser(prog='amptk-keep_samples.py',
                                     description='''Script parses AMPtk de-multiplexed FASTQ file and keeps those sequences with barocde names in list ''',
                                     epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)

    parser.add_argument('-i', '--input', required=True,
                        help='Input AMPtk demux FASTQ')
    parser.add_argument('-l', '--list', nargs='+',
                        help='Input list of (BC) names to keep')
    parser.add_argument('-t', '--threshold', type=int,
                        help='Keep samples with more reads than threshold')
    parser.add_argument('-f', '--file',
                        help='File containing list of names to keep')
    parser.add_argument('-o', '--out', required=True,
                        help='Output name')
    parser.add_argument('--format', default='fastq',
                        choices=['fastq', 'fasta'],
                        help='format of output file')
    args = parser.parse_args(args)

    keepers = []
    if args.threshold:
        print("Keeping samples with more than %i reads" % args.threshold)
        BC_counts = countBarcodes(args.input)
        for k, v in list(BC_counts.items()):
            if int(v) >= args.threshold:
                if not k in keepers:
                    keepers.append(k)

    if args.file:
        # load in list of sample names to keep
        with open(args.file, 'r') as input:
            lines = [line.rstrip('\n') for line in input]
        keepers += lines

    if args.list:
        lines = args.list
        keepers += lines

    # make sure it is a set, faster lookup
    keep_list = set(keepers)
    print("Keeping %i samples" % len(keep_list))

    # rename to base
    if args.out.endswith('.gz'):
        outfile = args.out.replace('.gz', '')
    else:
        outfile = args.out
    # run filtering
    keep_count, total_count = filter_sample(args.input, keep_list, outfile,
                                            format=args.format)
    # compress and clean
    if args.out.endswith('.gz'):  # compress in place
        amptklib.Fzip_inplace(outfile)

    print("Kept %i reads out of %i total reads" % (keep_count, total_count))


if __name__ == "__main__":
    main(sys.argv[1:])
