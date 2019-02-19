#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import re
import argparse
import logging
import subprocess
import inspect
import codecs
import unicodedata
import shutil
import multiprocessing
import datetime
from Bio import SeqIO
from amptk import amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

class col(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

latin_dict = {
u"¡": u"!", u"¢": u"c", u"£": u"L", u"¤": u"o", u"¥": u"Y",
u"¦": u"|", u"§": u"S", u"¨": u"`", u"©": u"c", u"ª": u"a",
u"«": u"<<", u"¬": u"-", u"­": u"-", u"®": u"R", u"¯": u"-",
u"°": u"o", u"±": u"+-", u"²": u"2", u"³": u"3", u"´": u"'",
u"µ": u"u", u"¶": u"P", u"·": u".", u"¸": u",", u"¹": u"1",
u"º": u"o", u"»": u">>", u"¼": u"1/4", u"½": u"1/2", u"¾": u"3/4",
u"¿": u"?", u"À": u"A", u"Á": u"A", u"Â": u"A", u"Ã": u"A",
u"Ä": u"A", u"Å": u"A", u"Æ": u"Ae", u"Ç": u"C", u"È": u"E",
u"É": u"E", u"Ê": u"E", u"Ë": u"E", u"Ì": u"I", u"Í": u"I",
u"Î": u"e", u"Ï": u"I", u"Ð": u"D", u"Ñ": u"N", u"Ò": u"O",
u"Ó": u"O", u"Ô": u"O", u"Õ": u"O", u"Ö": u"O", u"×": u"*",
u"Ø": u"O", u"Ù": u"U", u"Ú": u"U", u"Û": u"U", u"Ü": u"U",
u"Ý": u"Y", u"Þ": u"p", u"ß": u"b", u"à": u"a", u"á": u"a",
u"â": u"a", u"ã": u"a", u"ä": u"a", u"å": u"a", u"æ": u"ae",
u"ç": u"c", u"è": u"e", u"é": u"e", u"ê": u"e", u"ë": u"e",
u"ì": u"i", u"í": u"i", u"î": u"i", u"ï": u"i", u"ð": u"d",
u"ñ": u"n", u"ò": u"o", u"ó": u"o", u"ô": u"o", u"õ": u"o",
u"ö": u"o", u"÷": u"/", u"ø": u"o", u"ù": u"u", u"ú": u"u",
u"û": u"u", u"ü": u"u", u"ý": u"y", u"þ": u"p", u"ÿ": u"y",
u"’": u"'", u"×": u"x"}

def latin2ascii(error):
    return latin_dict[error.object[error.start]], error.start+1

#register codecs for latin conversion for strange characters
codecs.register_error('latin2ascii', latin2ascii)

def dereplicate(input, output, args=False):
    seqs = {}
    with open(output, 'w') as out:
        with open(input, 'r') as in_file:
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
                        if args.lca:
                            if newTaxLen == 1: #this means we have only a single level of taxonomy, so just move to next record
                                continue
                            if newTax[-1] == oldTax[-1]: #so taxonomy is the same, keep current value in dict and move to next record
                                continue
                            else: #loop backwards through tax string find last common ancestor
                                amptklib.log.debug("ERROR: %s and %s have identical sequences, but taxonomy doesn't agree" % (','.join(oldTax), ','.join(newTax)))
                                lca = 0
                                for num in range(1,newTaxLen+1):
                                    if newTax[-num] == oldTax[-num]:
                                        lca = num-1
                                        break
                                consensusTax = oldID+';tax='+','.join(oldTax[:-lca])
                                amptklib.log.debug("setting taxonomy to %s" % (consensusTax))
                                seqs[sequence] = consensusTax
        #now write to file     
        for key,value in seqs.items():
            out.write('>'+value+'\n'+key+'\n')
            
def uniteTax(taxString):
    '''
    function to parse UNITE taxonomy string, return taxstring
    '''
    tax, unite, reps, GenusSpecies, gbID = (None,)*5
    fields = taxString.split("|")
    for i in fields:
        if i.startswith("k__"):
            tax = i
        elif i.startswith("SH"):
            unite = i
        elif i.startswith("re"):
            reps = i
        elif '_' in i:
            GenusSpecies = i
        else:
            gbID = i
    if not unite:
        unite = 'NA'
    if not gbID:
       return None
    if tax:
        taxonomy = re.sub(";", ",", tax)
        taxonomy = re.sub("__", ":", taxonomy)
        tf = taxonomy.split(",")
        k = tf[0]
        k = re.sub('_', ' ', k)
        p = tf[1]
        p = re.sub('_', ' ', p)
        c = tf[2]
        c = re.sub('_', ' ', c)
        o = tf[3]
        o = re.sub('_', ' ', o)
        f = tf[4]
        f = re.sub('_', ' ', f)
        g = tf[5]
        g = re.sub('_', ' ', g)
        s = tf[6]
        s = re.sub('[(].*$','',s)
        s = re.sub('_', ' ', s)
        s = re.sub('\.', '', s)
        test_species = s.split(' ')
        if len(test_species) < 2:
            s = 's:'
        reformat_tax = []
        removal = ("unidentified", "Incertae", "uncultured", "Group", "incertae")
        sp_removal = (" sp", "_sp", "uncultured", "isolate", "mycorrhizae", "vouchered", "fungal", "basidiomycete", "ascomycete", "fungus", "symbiont")
        if not any(x in k for x in removal) and not amptklib.number_present(k):
            reformat_tax.append(k)
        if not any(x in p for x in removal) and not amptklib.number_present(p):
            reformat_tax.append(p)
        if not any(x in c for x in removal) and not amptklib.number_present(c):
            reformat_tax.append(c)
        if not any(x in o for x in removal) and not amptklib.number_present(o):
            reformat_tax.append(o)
        if not any(x in f for x in removal) and not amptklib.number_present(f):
            reformat_tax.append(f)
        if not any(x in g for x in removal) and not amptklib.number_present(g):
            reformat_tax.append(g)
        if not any(x in s for x in sp_removal):
            reformat_tax.append(s)
        taxReformatString = gbID+'|'+unite+";tax="+",".join(reformat_tax)
        taxReformatString = re.sub(",s:$", "", taxReformatString)
        taxReformatString = re.sub("=s:$", "=", taxReformatString)
        taxReformatString = taxReformatString.strip()
        return taxReformatString
    else:
        return None
    
def rdpTax(taxString):
    '''
    function to parse RDP taxonomy string, return seqrecord
    '''
    temp = taxString.split("\t")
    taxLevels = temp[-1]
    split_temp = temp[0].split(";")
    ID = split_temp[0].split(" ")[0]
    s = "s:" + split_temp[0].split(" ", 1)[-1]
    s = re.sub('[(].*$','',s)
    s = re.sub(',', '_', s)
    s = re.sub('\.', '', s)
    test_species = s.split(' ')
    if len(test_species) < 2:
        s = 's:'
    split_tax = taxLevels.split(";")
    if "domain" in split_tax:
        ki = split_tax.index("domain") -1
        k = "k:" + split_tax[ki]
        k = k.replace('"','')
        k = k.split(" ")[0]
    else:
        k = ""
    if "phylum" in split_tax:
        pi = split_tax.index("phylum") -1
        p = "p:" + split_tax[pi]
        p = p.replace('"','')
        p = p.split(" ")[0]
    else:
        p = ""
    if "class" in split_tax:
        ci = split_tax.index("class") -1
        c = "c:" + split_tax[ci]
        c = c.replace('"','')
        c = c.split(" ")[0]
    else:
        c = ""
    if "order" in split_tax:
        oi = split_tax.index("order") -1
        o = "o:" + split_tax[oi]
        o = o.replace('"','')
        o = o.split(" ")[0]
    else:
        o = ""
    if "family" in split_tax:
        fi = split_tax.index("family") -1
        f = "f:" + split_tax[fi]
        f = f.replace('"','')
        f = f.split(" ")[0]
    else:
        f = ""
    if "genus" in split_tax:
        gi = split_tax.index("genus") -1
        g = "g:" + split_tax[gi]
        g = g.replace('"','')
        g = g.split(" ")[0]
    else:
        g = ""
    reformat_tax = []
    removal = ("unidentified", "Incertae", "uncultured", "Group", "incertae", "Chloroplast", "unclassified", "Family", "leaf")
    sp_removal = (" sp", "_sp", "uncultured", "isolate", "mycorrhizae", "vouchered", "fungal", "basidiomycete", "ascomycete", "fungus", "symbiont", "unclassified", "unidentified", "bacterium", "phytoplasma", "genogroup")
    if not any(x in k for x in removal) and k != "" and not amptklib.number_present(k):
        reformat_tax.append(k)
    if not any(x in p for x in removal) and p != "" and not amptklib.number_present(p) and p != 'p:Fungi':
        reformat_tax.append(p)
    if not any(x in c for x in removal) and c != "" and not amptklib.number_present(c) and c != 'c:Fungi':
        reformat_tax.append(c)
    if not any(x in o for x in removal) and o != "" and not amptklib.number_present(o) and o != 'o:Fungi':
        reformat_tax.append(o)
    if not any(x in f for x in removal) and f != "" and not amptklib.number_present(f) and f != 'f:Fungi':
        reformat_tax.append(f)
    if not any(x in g for x in removal) and g != "" and not amptklib.number_present(g) and g != 'g:Fungi':
        reformat_tax.append(g)
    if not any(x in s for x in sp_removal):
        if ' cf ' in s or ' aff ' in s:
            chop = s.split(' ')
            if len(chop) > 3:
                s = chop[0]+ ' ' + chop[2]
            else:
                s = 's:'
        if ' I' in s:
            s = s.replace(" II", "")
            s = s.replace(" I", "")
        reformat_tax.append(s)
    taxReformatString = ID+";tax="+",".join(reformat_tax)
    taxReformatString = re.sub(",s:$", "", taxReformatString)
    taxReformatString = taxReformatString.strip()
    return taxReformatString        
      
def stripPrimer(input, args=False):
    base = os.path.basename(input).split('.')[0]
    StripOut = os.path.join(folder, base+'.extracted.fa')
    ErrorOut = os.path.join(folder, base+'.errors.fa')
    with open(StripOut, 'w') as outputfile:
        with open(ErrorOut, 'w') as errorfile:
            with open(input, 'r') as infile:
                for rec in SeqIO.parse(infile, 'fasta'):
                    orig_id = rec.description
                    if len(rec.seq) < args.min_len:
                        errorfile.write('>ERROR:LENGTH|%s\n%s\n' % (orig_id, str(rec.seq)))
                        continue
                    if args.utax == 'unite2utax':
                        newHeader = uniteTax(orig_id)
                    elif args.utax == 'rdp2utax':
                        newHeader = rdpTax(orig_id)
                    else:
                        newHeader = orig_id
                    #quality check the new header
                    if not newHeader:
                        errorfile.write('>ERROR:NO_ID|%s\n%s\n' % (orig_id, str(rec.seq)))
                        continue
                    if newHeader.strip().endswith('tax='):
                        errorfile.write('>ERROR:NO_TAX|%s\n%s\n' % (orig_id, str(rec.seq)))
                        continue
                    #rename header
                    rec.id = newHeader
                    rec.name = ''
                    rec.description = ''
                    #make sequence a string for processing
                    Seq = str(rec.seq).upper()
                    #remove any terminal N's
                    Seq = Seq.strip('N')
                    #make sure alignments are resetting
                    StripSeq, ForCutPos, RevCutPos, ForCutPos1, RevCutPos1 = (None,)*5
                    if not args.trimming:
                        #look for forward primer
                        ForCutPos = amptklib.findFwdPrimer(FwdPrimer, Seq, args.primer_mismatch, amptklib.degenNucSimple)
                        #align reverse primer, trim if found
                        RevCutPos = amptklib.findRevPrimer(RevPrimer, Seq, args.primer_mismatch, amptklib.degenNucSimple)
                        #now check if either were found
                        if ForCutPos and RevCutPos: #both were found, then trim and move on
                            StripSeq = Seq[ForCutPos:RevCutPos]
                        elif ForCutPos and not RevCutPos:
                            if args.keep == 'rev':
                                errorfile.write('>ERROR:NO_PRIMER_MATCH|%s\n%s\n' % (orig_id, str(rec.seq)))
                                continue
                            else:
                                StripSeq = Seq[ForCutPos:]
                        elif not ForCutPos and RevCutPos:
                            if args.keep == 'for':
                                errorfile.write('>ERROR:NO_PRIMER_MATCH|%s\n%s\n' % (orig_id, str(rec.seq)))
                                continue
                            else:
                                StripSeq = Seq[:RevCutPos]
                        elif not ForCutPos and not RevCutPos:
                            #neither is found, try reverse complementing
                            RevSeq = amptklib.RevComp(Seq)
                            ForCutPos1 = amptklib.findFwdPrimer(FwdPrimer, RevSeq, args.primer_mismatch, amptklib.degenNucSimple)
                            RevCutPos1 = amptklib.findRevPrimer(RevPrimer, RevSeq, args.primer_mismatch, amptklib.degenNucSimple)
                            #now check if either were found
                            if ForCutPos1 and RevCutPos1: #both were found, then trim and move on
                                StripSeq = RevSeq[ForCutPos1:RevCutPos1]
                            elif ForCutPos1 and not RevCutPos1:
                                if args.keep == 'rev':
                                    errorfile.write('>ERROR:NO_PRIMER_MATCH|%s\n%s\n' % (orig_id, str(rec.seq)))
                                    continue
                                else:
                                    StripSeq = RevSeq[ForCutPos1:]
                            elif not ForCutPos1 and RevCutPos1:
                                if args.keep == 'for':
                                    errorfile.write('>ERROR:NO_PRIMER_MATCH|%s\n%s\n' % (orig_id, str(rec.seq)))
                                    continue
                                else:
                                    StripSeq = RevSeq[:RevCutPos1]
                            elif not ForCutPos1 and not RevCutPos1: #neither is found, if keep all then return original sequence
                                if args.keep == 'none':
                                    StripSeq = Seq
                                else:
                                    errorfile.write('>ERROR:NO_PRIMER_MATCH|%s\n%s\n' % (orig_id, str(rec.seq)))
                                    continue 
                    else:
                        StripSeq = Seq                        
                    #check for ambiguous bases
                    if args.drop_ns != 0 and 'N'*args.drop_ns in StripSeq:
                        errorfile.write('>ERROR:AMBIGUOUS|%s\n%s\n' % (orig_id, Seq))
                        continue
                    #truncate
                    if args.trunclen:
                        StripSeq = StripSeq[:args.trunclen]
                    #check length
                    SeqLength = len(StripSeq)
                    if SeqLength >= args.min_len:
                        outputfile.write('>%s\n%s\n' % (rec.id, StripSeq))
                    else:
                        errorfile.write('>ERROR:LENGTH=%i|%s\n%s\n' % (SeqLength, orig_id, Seq))

def makeDB(input, args=False):
    basename = input.replace('.extracted.fa', '')
    db_details = basename + '.udb.txt'
    usearch_db = basename + '.udb'
    Total = amptklib.countfasta(input)
    if not ':' in args.source:
        source = args.source
        version = ''
    else:
        source, version = args.source.split(':')
    if args.trimming:
        args.F_primer = 'None'
        args.R_primer = 'None'
    if args.create_db == 'utax':
        databaseString = 'utax'
    elif args.create_db == 'usearch':
        databaseString = 'vsearch'
    elif args.create_db == 'sintax':
        databaseString = 'sintax'
    db_string = '{:} {:} {:} {:} {:} {:} {:} {:}'.format(databaseString, os.path.basename(args.fasta), args.F_primer, args.R_primer, str(Total), source, version, today)
    with open(db_details, 'w') as details:
        details.write(db_string)
    report = basename + '.report.txt'

    if args.create_db == 'utax':
        #create log file for this to troubleshoot
        utax_log = basename + '.utax.log'
        if os.path.isfile(utax_log):
            os.remove(utax_log)
        amptklib.log.info("Creating UTAX Database, this may take awhile")
        amptklib.log.debug("%s -makeudb_utax %s -output %s -report %s -utax_trainlevels kpcofgs -utax_splitlevels NVkpcofgs -notrunclabels" % (usearch, input, usearch_db, report))
        with open(utax_log, 'w') as utaxLog:
            subprocess.call([usearch, '-makeudb_utax', input, '-output', usearch_db, '-report', report, '-utax_trainlevels', args.utax_trainlevels, '-utax_splitlevels', args.utax_splitlevels, '-notrunclabels'], stdout = utaxLog, stderr = utaxLog)

        #check if file is actually there
        if os.path.isfile(usearch_db):
            amptklib.log.info("Database %s created successfully" % usearch_db)
        else:
            amptklib.log.error("There was a problem creating the DB, check the UTAX log file %s" % utax_log)

    elif args.create_db == 'usearch': #note we will use vsearch here
        #create log file for this to troubleshoot
        usearch_log = basename + '.usearch.log'
        if os.path.isfile(usearch_log):
            os.remove(usearch_log)
        amptklib.log.info("Creating USEARCH Database (VSEARCH)")
        amptklib.log.debug("vsearch --makeudb_usearch %s --output %s --notrunclabels" % (input, usearch_db))
        with open(usearch_log, 'w') as logfile:
            subprocess.call(['vsearch', '--makeudb_usearch', input, '--output', usearch_db, '--notrunclabels'], stdout = logfile, stderr = logfile)
        if os.path.isfile(usearch_db):
            amptklib.log.info("Database %s created successfully" % usearch_db)
        else:
            amptklib.log.error("There was a problem creating the DB, check the log file %s" % usearch_log)
    elif args.create_db == 'sintax':
        sintax_log = basename + '.sintax.log'
        if os.path.isfile(sintax_log):
            os.remove(sintax_log)
        amptklib.log.info("Creating SINTAX Database, this may take awhile")
        amptklib.log.debug("%s -makeudb_sintax %s -output %s" % (usearch, input, usearch_db))
        with open(sintax_log, 'w') as utaxLog:
            subprocess.call([usearch, '-makeudb_sintax', input, '-output', usearch_db, '-notrunclabels'], stdout = utaxLog, stderr = utaxLog)
        #check if file is actually there
        if os.path.isfile(usearch_db):
            amptklib.log.info("Database %s created successfully" % usearch_db)
        else:
            amptklib.log.error("There was a problem creating the DB, check the UTAX log file %s" % utax_log)
            
def decodeFasta(input, output):
    '''
    the taxonomy strings from UNITE database contain invalid characters, try to fix them
    '''
    with open(output, 'w') as outfile:
        with open(input, 'r') as infile:
            for rec in SeqIO.parse(infile, 'fasta'):
                try:
                    tmpID = rec.description.decode('utf-8')
                except AttributeError:
                    tmpID = rec.description
                cleanID = tmpID.encode('ascii', 'latin2ascii').decode()
                outfile.write('>{:}\n{:}\n'.format(cleanID, str(rec.seq)))

def main(args):
    global FwdPrimer, RevPrimer, today, folder, usearch
    parser=argparse.ArgumentParser(prog='amptk-extract_region.py', usage="%(prog)s [options] -f <FASTA File>",
        description='''Script searches for primers and removes them if found.  Useful for trimming a reference dataset for assigning taxonomy after OTU clustering.  It is also capable of reformatting UNITE taxonomy fasta headers to be compatible with UTAX and creating USEARCH/UTAX UBD databases for assigning taxonomy.''',
        epilog="""Written by Jon Palmer (2015-2017) nextgenusfs@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('-i','--fasta', dest='fasta', required=True, help='FASTA input')
    parser.add_argument('-o','--out', dest='out', help='Base Name Output files')
    parser.add_argument('-f','--fwd_primer', dest='F_primer', default='fITS7', help='Forward primer (fITS7)')
    parser.add_argument('-r','--rev_primer', dest='R_primer', default='ITS4', help='Reverse primer (ITS4)')
    parser.add_argument('--skip_trimming', dest='trimming', action='store_true', help='Skip Primer trimming (not recommended)')
    parser.add_argument('--format', dest='utax', default='unite2utax', choices=['unite2utax', 'rdp2utax', 'off'], help='Reformat FASTA headers for UTAX')
    parser.add_argument('--lca', action='store_true', help='Run LCA (last common ancestor) for dereplicating taxonomy')
    parser.add_argument('--trunclen', type=int, help='Truncate reads to length')
    parser.add_argument('--subsample', type=int, help='Random subsample')
    parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
    parser.add_argument('--max_len', default=1200, type=int, help='Maximum read length to keep')
    parser.add_argument('--drop_ns', dest='drop_ns', type=int, default=8, help="Drop Seqeunces with more than X # of N's")
    parser.add_argument('--create_db', dest='create_db', choices=['utax', 'usearch', 'sintax'], help="Create USEARCH DB")
    parser.add_argument('--primer_required', dest='keep', default='for', choices=['none', 'for', 'rev'], help="Keep Seq if primer found Default: for")
    parser.add_argument('--derep_fulllength', action='store_true', help="De-replicate sequences. Default: off")
    parser.add_argument('--primer_mismatch', default=2, type=int, help="Max Primer Mismatch")
    parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: auto")
    parser.add_argument('--install', action='store_true', help="Install into AMPtk database")
    parser.add_argument('--source', default=':', help="DB source and version separated by :")
    parser.add_argument('--utax_trainlevels', default='kpcofgs', help="UTAX training parameters")
    parser.add_argument('--utax_splitlevels', default='NVkpcofgs', help="UTAX training parameters")
    parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
    parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
    args=parser.parse_args(args)
    
    parentdir = os.path.join(os.path.dirname(amptklib.__file__))

    #get basename if not args.out passed
    if args.out:
        base = args.out
    else:
        if 'otus' in args.input:
            base = os.path.basename(args.fasta).split('.otus')[0]
        else:
            base = os.path.basename(args.fasta).split('.f')[0]
    
    if args.install: #set path to the AMPtk database location
        base = os.path.join(parentdir, 'DB', base)
        
    #remove logfile if exists
    log_name = base + '.log'
    if os.path.isfile(log_name):
        os.remove(log_name)

    amptklib.setupLogging(log_name)
    FNULL = open(os.devnull, 'w')
    cmd_args = " ".join(sys.argv)+'\n'
    amptklib.log.debug(cmd_args)
    print("-------------------------------------------------------")
    amptklib.SystemInfo()

    today = datetime.datetime.today().strftime('%Y-%m-%d')

    #Do a version check
    usearch = args.usearch
    amptklib.versionDependencyChecks(usearch)
    
    amptklib.log.info('Base name set to: {:}'.format(base))

    #look up primer db otherwise default to entry
    if args.F_primer in amptklib.primer_db:
        FwdPrimer = amptklib.primer_db.get(args.F_primer)
    else:
        FwdPrimer = args.F_primer
    if args.R_primer in amptklib.primer_db:
        RevPrimer = amptklib.primer_db.get(args.R_primer)
    else:
        RevPrimer = args.R_primer

    #reverse complement the reverse primer and get forward length for trim function
    RevPrimer = amptklib.RevComp(RevPrimer)
    fwdLen = len(FwdPrimer)

    if not args.trimming:
        amptklib.log.info("Searching for primers, this may take awhile: Fwd: %s  Rev: %s" % (FwdPrimer, RevPrimer))
    else:
        amptklib.log.info("Working on file: %s" % args.fasta)

    #get number of cpus
    if args.cpus:
        cpus = args.cpus
    else:
        cpus = amptklib.getCPUS()

    #create temp directory
    pid = os.getpid()
    folder = 'amptk_tmp_' + str(pid)
    if not os.path.exists(folder):
        os.makedirs(folder)

    decodedFasta = os.path.join(folder, 'cleaned.fa')
    decodeFasta(args.fasta, decodedFasta)
    SeqCount = amptklib.countfasta(decodedFasta)
    amptklib.log.info('{0:,}'.format(SeqCount) + ' records loaded')
    #if only 1 cpu just process data
    if cpus == 1:
        stripPrimer(decodedFasta, args=args)
    else:
        amptklib.log.info("Using %i cpus to process data" % cpus)

        #now split it into chunks (as many cpus as are queried)
        amptklib.split_fasta(decodedFasta, folder, cpus*2)

        #get list of files
        file_list = []
        for file in os.listdir(folder):
            if file.endswith(".fasta"):
                file = os.path.join(folder, file)
                file_list.append(file)

        #finally process reads over number of cpus
        amptklib.runMultiProgress(stripPrimer, file_list, cpus, args=args)

    #now concatenate outputs together
    OutName = base + '.extracted.fa'
    ErrorName = base + '.errors.fa'
    with open(OutName, 'w') as outfile:
        with open(ErrorName, 'w') as outfile2:
            for filename in os.listdir(os.path.join(folder)):
                if filename.endswith('.extracted.fa'):
                    if filename == OutName:
                        continue
                    with open(os.path.join(folder, filename), 'r') as readfile:
                        shutil.copyfileobj(readfile, outfile)
                if filename.endswith('.errors.fa'):
                    if filename == ErrorName:
                        continue
                    with open(os.path.join(folder, filename), 'r') as readfile:
                        shutil.copyfileobj(readfile, outfile2)

    if not args.debug:
        #clean up tmp folder
        shutil.rmtree(folder)

    #parse stats from the error fasta file
    noID = 0
    ambig = 0
    tooShort = 0
    noPrimer = 0
    noTax = 0
    with open(ErrorName, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                Err = line.split('|')[0]
                Err = Err.replace('>ERROR:', '')
                if Err == 'NO_ID':
                    noID += 1
                elif Err == 'AMBIGUOUS':
                    ambig += 1
                elif 'LENGTH' in Err:
                    tooShort += 1
                elif Err == 'NO_PRIMER_MATCH':
                    noPrimer += 1
                elif Err == 'NO_TAX':
                    noTax += 1

    Passed = amptklib.countfasta(os.path.abspath(OutName))
    PassedPct = Passed / SeqCount
    amptklib.log.info('{:,} records passed ({:.2%})'.format(Passed, PassedPct))
    amptklib.log.info('Errors: {:,} no taxonomy info, {:,} no ID, {:,} length out of range, {:,} too many ambiguous bases, {:,} no primers found'.format(
            noTax, noID, tooShort, ambig, noPrimer))

    #dereplicate if argument passed
    if args.derep_fulllength:
        amptklib.log.info("Now dereplicating sequences (collapsing identical sequences)")
        derep_tmp = base + '.derep.extracted.fa'
        os.rename(OutName, derep_tmp)
        dereplicate(derep_tmp, OutName, args=args)
        Total = amptklib.countfasta(OutName)
        amptklib.log.info('{0:,}'.format(Total) + ' records passed (%.2f%%)' % (Total*100.0/SeqCount))
        os.remove(derep_tmp)

    #subsample if argument passed
    if args.subsample:
        amptklib.log.info("Random subsampling to %i records" % args.subsample)
        subsample = base + '.subsample.tmp'
        os.rename(OutName, subsample)
        cmd = ['vsearch', '--fastx_subsample', subsample, '--sample_size', str(args.subsample), '--fastaout', OutName, '--notrunclabels']
        amptklib.runSubprocess(cmd, amptklib.log)
        os.remove(subsample)

    #create DB if argument passed  
    if args.create_db:
        makeDB(OutName, args=args)
    print("-------------------------------------------------------")
            
if __name__ == "__main__":
    main(args)
