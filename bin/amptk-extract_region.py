#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, re, argparse, logging, subprocess, inspect, codecs, unicodedata, shutil, glob, multiprocessing
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.amptklib as amptklib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='amptk-extract_region.py', usage="%(prog)s [options] -f <FASTA File>",
    description='''Script searches for primers and removes them if found.  Useful for trimming a reference dataset for assigning taxonomy after OTU clustering.  It is also capable of reformatting UNITE taxonomy fasta headers to be compatible with UTAX and creating USEARCH/UTAX UBD databases for assigning taxonomy.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fasta', dest='fasta', required=True, help='FASTA input')
parser.add_argument('-o','--out', dest='out', default='amptk', help='Base Name Output files')
parser.add_argument('-f','--fwd_primer', dest='F_primer', default='fITS7', help='Forward primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest='R_primer', default='ITS4', help='Reverse primer (ITS4)')
parser.add_argument('--skip_trimming', dest='trimming', action='store_true', help='Skip Primer trimming (not recommended)')
parser.add_argument('--format', dest='utax', default='unite2utax', choices=['unite2utax', 'rdp2utax', 'off'], help='Reformat FASTA headers for UTAX')
parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
parser.add_argument('--max_len', default=1200, type=int, help='Maximum read length to keep')
parser.add_argument('--drop_ns', dest='drop_ns', type=int, default=8, help="Drop Seqeunces with more than X # of N's")
parser.add_argument('--create_db', dest='create_db', choices=['utax', 'usearch'], help="Create USEARCH DB")
parser.add_argument('--keep_all', dest='keep_all', action='store_true', help="Keep Seq if For primer not found Default: off")
parser.add_argument('--derep_fulllength', action='store_true', help="De-replicate sequences. Default: off")
parser.add_argument('--primer_mismatch', default=4, help="Max Primer Mismatch")
parser.add_argument('--cpus', type=int, help="Number of CPUs. Default: max")
parser.add_argument('--utax_trainlevels', default='kpcofgs', help="UTAX training parameters")
parser.add_argument('--utax_splitlevels', default='NVkpcofgs', help="UTAX training parameters")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH8 EXE')
args=parser.parse_args()

#look up primer db otherwise default to entry
primer_db = {'fITS7': 'GTGARTCATCGAATCTTTG', 'ITS4': 'TCCTCCGCTTATTGATATGC', 'ITS1-F': 'CTTGGTCATTTAGAGGAAGTAA', 'ITS2': 'GCTGCGTTCTTCATCGATGC', 'ITS3': 'GCATCGATGAAGAACGCAGC', 'ITS4-B': 'CAGGAGACTTGTACACGGTCCAG', 'ITS1': 'TCCGTAGGTGAACCTGCGG', 'LR0R': 'ACCCGCTGAACTTAAGC', 'LR2R': 'AAGAACTTTGAAAAGAG', 'JH-LS-369rc': 'CTTCCCTTTCAACAATTTCAC', '16S_V3': 'CCTACGGGNGGCWGCAG', '16S_V4': 'GACTACHVGGGTATCTAATCC'}
if args.F_primer in primer_db:
    FwdPrimer = primer_db.get(args.F_primer)
else:
    FwdPrimer = args.F_primer
if args.R_primer in primer_db:
    RevPrimer = primer_db.get(args.R_primer)
else:
    RevPrimer = args.R_primer

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
u"’":u"'", u"×":u"x"}

def latin2ascii(error):
    return latin_dict[error.object[error.start]], error.start+1
codecs.register_error('latin2ascii', latin2ascii)

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
                    newTax = rec.description.split(',')
                    oldTax = seqs.get(sequence).split(',')
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
                            amptklib.log.debug("ERROR: %s and %s have identical sequences, but taxonomy doesn't agree" % (','.join(oldTax), ','.join(newTax)))
                            lca = 0
                            for num in range(1,newTaxLen+1):
                                if newTax[-num] == oldTax[-num]:
                                    lca = num-1
                                    break
                            consensusTax = ','.join(oldTax[:-lca])
                            amptklib.log.debug("setting taxonomy to %s" % (consensusTax))
                            seqs[sequence] = consensusTax

        #now write to file     
        for key,value in seqs.iteritems():
            out.write('>'+value+'\n'+key+'\n')

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def stripPrimer(records):
    for rec in records:
        if args.utax == 'unite2utax':
            latin = unicode(rec.description, 'utf-8')
            test = latin.encode('ascii', 'latin2ascii')
            fields = test.split("|")
            for i in fields:
                if i.startswith("k__"):
                    tax = i
                elif i.startswith("SH"):
                    unite = i
                elif i.startswith("re"):
                    reps = i
                else:
                    gbID = i
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
            if not any(x in k for x in removal):
                reformat_tax.append(k)
            if not any(x in p for x in removal):
                reformat_tax.append(p)
            if not any(x in c for x in removal):
                reformat_tax.append(c)
            if not any(x in o for x in removal):
                reformat_tax.append(o)
            if not any(x in f for x in removal):
                reformat_tax.append(f)
            if not any(x in g for x in removal):
                reformat_tax.append(g)
            if not any(x in s for x in sp_removal):
                reformat_tax.append(s)
            rec.id = gbID+";tax="+",".join(reformat_tax)
            rec.id = re.sub(",s:$", "", rec.id)
            rec.id = re.sub("=s:$", "=", rec.id)
            if rec.id.endswith(";tax="): #if there is no taxonomy, get rid of it
                rec.id = ""
            rec.name = ""
            rec.description = ""
        elif args.utax == 'rdp2utax':
            latin = unicode(rec.description, 'utf-8')
            test = latin.encode('ascii', 'latin2ascii')
            temp = test.split("\t")
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
            removal = ("unidentified", "Incertae", "uncultured", "Group", "incertae", "Chloroplast", "unclassified", "Family")
            sp_removal = (" sp", "_sp", "uncultured", "isolate", "mycorrhizae", "vouchered", "fungal", "basidiomycete", "ascomycete", "fungus", "symbiont", "unclassified", "unidentified", "bacterium", "phytoplasma")
            if not any(x in k for x in removal) and k != "":
                reformat_tax.append(k)
            if not any(x in p for x in removal) and p != "":
                reformat_tax.append(p)
            if not any(x in c for x in removal) and c != "":
                reformat_tax.append(c)
            if not any(x in o for x in removal) and o != "":
                reformat_tax.append(o)
            if not any(x in f for x in removal) and f != "":
                reformat_tax.append(f)
            if not any(x in g for x in removal) and g != "":
                reformat_tax.append(g)
            if not any(x in s for x in sp_removal):
                reformat_tax.append(s)
            rec.id = ID+";tax="+",".join(reformat_tax)
            rec.id = re.sub(",s:$", "", rec.id)
            if rec.id.endswith(";tax="): #if there is no taxonomy, get rid of it
                rec.id = ""
            rec.name = ""
            rec.description = ""
        if not args.trimming:
            Seq = rec.seq
            MAX_PRIMER_MISMATCHES = int(args.primer_mismatch)
            revPrimer = revcomp_lib.RevComp(RevPrimer)
            BestPosFor, BestDiffsFor = primer.BestMatch2(Seq, FwdPrimer, MAX_PRIMER_MISMATCHES)
            if BestDiffsFor < MAX_PRIMER_MISMATCHES:
                if BestPosFor > 0:
                    stripfwdlen = fwdLen + BestPosFor
                    StripSeq = Seq[stripfwdlen:]

                    #now look for reverse
                    BestPosRev, BestDiffsRev = primer.BestMatch2(StripSeq, revPrimer, MAX_PRIMER_MISMATCHES)
                    if BestDiffsRev < MAX_PRIMER_MISMATCHES:
                        StrippedSeq = StripSeq[:BestPosRev]
                    else:
                        StrippedSeq = StripSeq
                    #after stripping primers, check for ambig bases
                    if args.drop_ns != 0 and 'N'*args.drop_ns in StrippedSeq:
                        continue
                    rec.seq = StrippedSeq
                    if rec.id != "" and rec.seq != "" and len(rec.seq) >= args.min_len and len(rec.seq) <= args.max_len:
                        yield rec
            else: #if can't find forward primer, try to reverse complement and look again
                RevSeq = revcomp_lib.RevComp(Seq)
                BestPosFor, BestDiffsFor = primer.BestMatch2(RevSeq, FwdPrimer, MAX_PRIMER_MISMATCHES)
                if BestDiffsFor < MAX_PRIMER_MISMATCHES:
                    if BestPosFor > 0:
                        stripfwdlen = fwdLen + BestPosFor
                        StripSeq = Seq[stripfwdlen:]

                        #now look for reverse
                        BestPosRev, BestDiffsRev = primer.BestMatch2(StripSeq, revPrimer, MAX_PRIMER_MISMATCHES)
                        if BestDiffsRev < MAX_PRIMER_MISMATCHES:
                            StrippedSeq = StripSeq[:BestPosRev]
                        else:
                            StrippedSeq = StripSeq
                        #after stripping primers, check for ambig bases
                        if args.drop_ns != 0 and 'N'*args.drop_ns in StrippedSeq:
                            continue
                        rec.seq = StrippedSeq
                        if rec.id != "" and rec.seq != "" and len(rec.seq) >= args.min_len and len(rec.seq) <= args.max_len:
                            yield rec
                else:
                    if args.keep_all:
                        StripSeq = Seq
                        #now look for reverse
                        BestPosRev, BestDiffsRev = primer.BestMatch2(StripSeq, revPrimer, MAX_PRIMER_MISMATCHES)
                        if BestDiffsRev < MAX_PRIMER_MISMATCHES:
                            StrippedSeq = StripSeq[:BestPosRev]
                        else:
                            StrippedSeq = StripSeq
                        #after stripping primers, check for ambig bases
                        if args.drop_ns != 0 and 'N'*args.drop_ns in StrippedSeq:
                            continue
                        rec.seq = StrippedSeq
                        if rec.id != "" and rec.seq != "" and len(rec.seq) >= args.min_len and len(rec.seq) <= args.max_len:
                            yield rec
        else:
            #check for ambig bases
            Seq = str(rec.seq)
            if args.drop_ns != 0 and 'N'*args.drop_ns in Seq:
                continue
            if rec.id != "" and rec.seq != "" and len(rec.seq) >= args.min_len and len(rec.seq) <= args.max_len:
                yield rec

def makeDB(input):
    #need usearch for this, test to make sure version is ok with utax
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        amptklib.log.error("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    version = usearch_test.split(" v")[1]
    majorV = version.split(".")[0]
    minorV = version.split(".")[1]
    if int(majorV) < 8 or (int(majorV) >= 8 and int(minorV) < 1):
        amptklib.log.warning("USEARCH version: %s detected you need v8.1.1756 or above" % usearch_test)
        os._exit(1)
    else:
        amptklib.log.info("USEARCH version: %s" % usearch_test)

    db_details = args.out + '.udb.txt'
    usearch_db = args.out + '.udb'
    if args.trimming:
        args.F_primer = 'None'
        args.R_primer = 'None'
    db_string = args.create_db + ' ' + args.fasta + ' ' + args.F_primer + ' ' + args.R_primer + ' ' + str(Total)
    with open(db_details, 'w') as details:
        details.write(db_string)
    report = args.out + '.report.txt'


    if args.create_db == 'utax':
        #create log file for this to troubleshoot
        utax_log = args.out + '.utax.log'
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

    if args.create_db == 'usearch':
        #create log file for this to troubleshoot
        usearch_log = args.out + '.usearch.log'
        if os.path.isfile(usearch_log):
            os.remove(usearch_log)
        amptklib.log.info("Creating USEARCH Database")
        amptklib.log.debug("%s -makeudb_usearch %s -output %s -notrunclabels" % (usearch, input, usearch_db))
        with open(usearch_log, 'w') as logfile:
            subprocess.call([usearch, '-makeudb_usearch', input, '-output', usearch_db, '-notrunclabels'], stdout = logfile, stderr = logfile)
        if os.path.isfile(usearch_db):
            amptklib.log.info("Database %s created successfully" % usearch_db)
        else:
            amptklib.log.error("There was a problem creating the DB, check the log file %s" % utax_log)

def batch_iterator(iterator, batch_size):
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def worker(input):
    output = input.split(".",-1)[0] + '.extracted.fa'
    with open(output, 'w') as o:
        with open(input, 'rU') as i:
            SeqRecords = SeqIO.parse(i, 'fasta')
            SeqIO.write(stripPrimer(SeqRecords), o, 'fasta')


#remove logfile if exists
log_name = args.out + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)

amptklib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
amptklib.log.debug(cmd_args)
print "-------------------------------------------------------"
amptklib.SystemInfo()

FileName = args.fasta
fwdLen = len(FwdPrimer)

if not args.trimming:
    amptklib.log.info("Searching for primers, this may take awhile: Fwd: %s  Rev: %s" % (args.F_primer, args.R_primer))
else:
    amptklib.log.info("Working on file: %s" % FileName)

if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

with open(FileName, 'rU') as input:
    SeqCount = countfasta(FileName)
    amptklib.log.info('{0:,}'.format(SeqCount) + ' records loaded')
    SeqRecords = SeqIO.parse(FileName, 'fasta')
    chunks = SeqCount / cpus + 1
    amptklib.log.info("Using %i cpus to process data" % cpus)
    #divide into chunks, store in tmp file
    pid = os.getpid()
    folder = 'amptk_tmp_' + str(pid)
    if not os.path.exists(folder):
        os.makedirs(folder)
    for i, batch in enumerate(batch_iterator(SeqRecords, chunks)) :
        filename = "chunk_%i.fa" % (i+1)
        tmpout = os.path.join(folder, filename)
        handle = open(tmpout, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()

file_list = []
for file in os.listdir(folder):
    if file.endswith(".fa"):
        file = os.path.join(folder, file)
        file_list.append(file)

p = multiprocessing.Pool(cpus)
for f in file_list:
    p.apply_async(worker, [f])
p.close()
p.join()

#now concatenate outputs together
OutName = args.out + '.extracted.fa'
with open(OutName, 'w') as outfile:
    for filename in glob.glob(os.path.join(folder,'*.extracted.fa')):
        if filename == OutName:
            continue
        with open(filename, 'rU') as readfile:
            shutil.copyfileobj(readfile, outfile)
#clean up tmp folder
shutil.rmtree(folder)

if args.derep_fulllength:
    Passed = countfasta(OutName)
    amptklib.log.info('{0:,}'.format(Passed) + ' records passed (%.2f%%)' % (Passed*100.0/SeqCount))
    amptklib.log.info("Now dereplicating sequences (remove if sequence and header identical)")
    derep_tmp = args.out + '.derep.extracted.fa'
    os.rename(OutName, derep_tmp)
    dereplicate(derep_tmp, OutName)
    Total = countfasta(OutName)
    amptklib.log.info('{0:,}'.format(Total) + ' records passed (%.2f%%)' % (Total*100.0/SeqCount))
    os.remove(derep_tmp)
else:
    Total = countfasta(OutName)
    amptklib.log.info('{0:,}'.format(Total) + ' records passed (%.2f%%)' % (Total*100.0/SeqCount))

if args.create_db:
    makeDB(OutName)

print "-------------------------------------------------------"