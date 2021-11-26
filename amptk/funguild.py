#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import argparse
import json
import os
import sys
import amptk.amptklib as lib
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

# simple script to replicate funguild methods which do not parse
# AMPtk taxonomy properly

def load_database(url):
    # will return dictionary from funguild
    db_content = json.load(urlopen(url))
    search_dict = {}
    for rec in db_content:
        search_dict[rec['taxon']] = rec
    return search_dict


def run_funguild(otutable, guildedtable,
                 dburl='https://mycoportal.org/fdex/services/api/db_return.php?dbReturn=Yes&pp=1'):
    # wrapper for running amptk funguild
    guild_db = load_database(dburl)
    # we can now open tax table, parse taxonomy, search agains guilds, and append in one go
    guild_headers = ['taxon', 'taxonomicLevel', 'trophicMode', 'guild',
                     'trait', 'growthForm', 'confidenceRanking', 'notes',
                     'citationSource']
    stats = {'dbsize': len(guild_db),
             'trophic-levels': len(set([v['trophicMode'] for k,v in guild_db.items()])),
             'guilds': len(set([v['guild'] for k,v in guild_db.items()])),
             'classified': {
                'otus': 0,
                'u': 0,
                'levels': {}
                }
             }
    with open(guildedtable, 'w') as outfile:
        with open(otutable, 'r') as infile:
            for line in infile:
                found = False
                line = line.rstrip()
                if line.startswith('\n'):
                    continue
                cols = line.split('\t')
                if line.startswith('#OTU'):
                    newline = cols + guild_headers
                    outfile.write('{}\n'.format('\t'.join(newline)))
                else:
                    stats['classified']['otus'] += 1
                    # taxonomy is in last column
                    fullTax = cols[-1].rsplit(';',1)[-1]
                    if fullTax in ['', 'No hit']:
                        newline = cols + ['',]*8
                    else:
                        taxList = fullTax.split(',')
                        # can now loop through the taxonomy list in reverse, stopping if find something
                        for tax in taxList[::-1]:
                            if found:
                                continue
                            level, taxname = tax.split(':')
                            if taxname in guild_db:
                                hit = guild_db[taxname]
                                newline = cols + [hit['taxon'], hit['taxonomicLevel'], hit['trophicMode'],
                                                hit['guild'], hit['growthForm'], hit['confidenceRanking'],
                                                hit['notes'], hit['citationSource']
                                                ]
                                if not level in stats['classified']:
                                    stats['classified'][level] = 1
                                else:
                                    stats['classified'][level] += 1
                                if not hit['trophicMode'] in stats['classified']['levels']:
                                    stats['classified']['levels'][hit['trophicMode']] = 1
                                else:
                                    stats['classified']['levels'][hit['trophicMode']] += 1
                                found = True
                    if not found:
                        stats['classified']['u'] += 1
                        newline = cols + ['',]*8
                    outfile.write('{}\n'.format('\t'.join(newline)))
    return stats


def main(args):
    parser=argparse.ArgumentParser(prog='amptk-funguild.py',
        description='''run FUNGuild analysis of OTU table''',
        epilog="""Written by Jon Palmer (2021) nextgenusfs@gmail.com""")
    parser.add_argument('-i','--input', required=True, help='OTU table with Taxonomy from AMPtk')
    parser.add_argument('-o','--out', help='Output OTU table')
    parser.add_argument('-d','--database', default='fungi', help='Option not used, here for compatibility')
    parser.add_argument('-u','--url',
                        default='https://mycoportal.org/fdex/services/api/db_return.php?dbReturn=Yes&pp=1',
                        help='FUNGuild API url')
    args=parser.parse_args(args)

    # initialize script and log
    if not args.out:
        args.out = '{}.funguild.txt'.format(args.input.rsplit('.', 1)[0])

    #remove logfile if exists
    log_name = '{}.log'.format(args.out.rsplit('.', 1)[0])
    if os.path.isfile(log_name):
        os.remove(log_name)

    lib.setupLogging(log_name)
    cmd_args = '{}\n'.format(" ".join(sys.argv))
    lib.log.debug(cmd_args)
    sys.stderr.write("-------------------------------------------------------\n")
    # check for taxonomy quickly
    lib.SystemInfo()
    with open(args.input) as f:
        first_line = f.readline()
        if not 'Taxonomy' in first_line:
            lib.log.error('Taxonomy field is not found in OTU table, exiting')
            sys.exit(1)

    # if above passes, then run the method
    lib.log.info('Downloading/parsing FUNGuild database from: {}'.format(args.url))
    stats = run_funguild(args.input, args.out, dburl=args.url)
    lib.log.info('Assigning functional guilds completed')
    lib.log.debug(stats)
    lib.log.info('FUNGuild databases consists of {:,} records, {:,} trophic-levels, {:,} guilds'.format(stats['dbsize'], stats['trophic-levels'], stats['guilds']))
    lib.log.info('Trophic-level assignment statistics for {:,} OTUS'.format(stats['classified']['otus']))
    for k,v in sorted(stats['classified']['levels'].items()):
        if not '-' in k:
            sys.stderr.write('{}\t\t{}\n'.format(k, v))
        else:
            sys.stderr.write('{}\t{}\n'.format(k, v))
    sys.stderr.write('Unclassified\t\t{}\n'.format(stats['classified']['u']))
    sys.stderr.write("-------------------------------------------------------\n")


if __name__ == "__main__":
    main(sys.argv[1:])