#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, re, argparse, logging, subprocess, inspect, codecs, unicodedata, shutil, multiprocessing
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
    epilog="""Written by Jon Palmer (2015-2017) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fasta', dest='fasta', required=True, help='FASTA input')
parser.add_argument('-o','--out', dest='out', default='amptk', help='Base Name Output files')
parser.add_argument('-f','--fwd_primer', dest='F_primer', default='fITS7', help='Forward primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest='R_primer', default='ITS4', help='Reverse primer (ITS4)')
parser.add_argument('--skip_trimming', dest='trimming', action='store_true', help='Skip Primer trimming (not recommended)')
parser.add_argument('--format', dest='utax', default='unite2utax', choices=['unite2utax', 'rdp2utax', 'off'], help='Reformat FASTA headers for UTAX')
parser.add_argument('--lca', action='store_true', help='Run LCA (last common ancestor) for dereplicating taxonomy')
parser.add_argument('--min_len', default=100, type=int, help='Minimum read length to keep')
parser.add_argument('--max_len', default=1200, type=int, help='Maximum read length to keep')
parser.add_argument('--drop_ns', dest='drop_ns', type=int, default=8, help="Drop Seqeunces with more than X # of N's")
parser.add_argument('--create_db', dest='create_db', choices=['utax', 'usearch'], help="Create USEARCH DB")
parser.add_argument('--keep_all', dest='keep_all', action='store_true', help="Keep Seq if For primer not found Default: off")
parser.add_argument('--derep_fulllength', action='store_true', help="De-replicate sequences. Default: off")
parser.add_argument('--primer_mismatch', default=2, help="Max Primer Mismatch")
parser.add_argument('--cpus', type=int, default=1, help="Number of CPUs. Default: 1")
parser.add_argument('--utax_trainlevels', default='kpcofgs', help="UTAX training parameters")
parser.add_argument('--utax_splitlevels', default='NVkpcofgs', help="UTAX training parameters")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
args=parser.parse_args()

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

#register codecs for latin conversion for strange characters
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
        for key,value in seqs.iteritems():
            out.write('>'+value+'\n'+key+'\n')
        
def stripPrimer(input):
    base = os.path.basename(input).split('.')[0]
    StripOut = os.path.join(folder, base+'.extracted.fa')
    ErrorOut = os.path.join(folder, base+'.errors.fa')
    with open(StripOut, 'w') as outputfile:
        with open(ErrorOut, 'w') as errorfile:
            with open(input, 'rU') as infile:
                for rec in SeqIO.parse(infile, 'fasta'):
                    orig_id = rec.description
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
                        if not any(x in k for x in removal) and k != "" and not amptklib.number_present(k):
                            reformat_tax.append(k)
                        if not any(x in p for x in removal) and p != "" and not amptklib.number_present(p):
                            reformat_tax.append(p)
                        if not any(x in c for x in removal) and c != "" and not amptklib.number_present(c):
                            reformat_tax.append(c)
                        if not any(x in o for x in removal) and o != "" and not amptklib.number_present(o):
                            reformat_tax.append(o)
                        if not any(x in f for x in removal) and f != "" and not amptklib.number_present(f):
                            reformat_tax.append(f)
                        if not any(x in g for x in removal) and g != "" and not amptklib.number_present(g):
                            reformat_tax.append(g)
                        if not any(x in s for x in sp_removal):
                            reformat_tax.append(s)
                        rec.id = ID+";tax="+",".join(reformat_tax)
                        rec.id = re.sub(",s:$", "", rec.id)
                        if rec.id.endswith(";tax="): #if there is no taxonomy, get rid of it
                            rec.id = ""
                        rec.name = ""
                        rec.description = ""                
                    if rec.id == "":
                        errorfile.write('>ERROR:NO_ID|%s\n%s\n' % (orig_id, str(rec.seq)))
                        continue
                    #make sequence a string for processing
                    Seq = str(rec.seq)
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
                            StripSeq = Seq[ForCutPos:]
                        elif not ForCutPos and RevCutPos:
                            StripSeq = Seq[:RevCutPos]
                        elif not ForCutPos and not RevCutPos:
                            #neither is found, try reverse complementing
                            RevSeq = revcomp_lib.RevComp(Seq)
                            ForCutPos1 = amptklib.findFwdPrimer(FwdPrimer, RevSeq, args.primer_mismatch, amptklib.degenNucSimple)
                            RevCutPos1 = amptklib.findRevPrimer(RevPrimer, RevSeq, args.primer_mismatch, amptklib.degenNucSimple)
                            #now check if either were found
                            if ForCutPos1 and RevCutPos1: #both were found, then trim and move on
                                StripSeq = RevSeq[ForCutPos1:RevCutPos1]
                            elif ForCutPos1 and not RevCutPos1:
                                StripSeq = RevSeq[ForCutPos1:]
                            elif not ForCutPos1 and RevCutPos1:
                                StripSeq = RevSeq[:RevCutPos1]
                            elif not ForCutPos1 and not RevCutPos1: #neither is found, if keep all then return original sequence
                                if args.keep_all:
                                    StripSeq = Seq
                                else:
                                    errorfile.write('>ERROR:NO_PRIMER_MATCH|%s\n%s\n' % (orig_id, str(rec.seq)))
                                    continue 
                    else:
                        StripSeq = Seq                             
                    #check for ambiguous bases
                    if args.drop_ns != 0 and 'N'*args.drop_ns in StripSeq:
                        errorfile.write('>ERROR:AMBIGUOUS|%s\n%s\n' % (orig_id, str(rec.seq)))
                        continue
                    #check length
                    SeqLength = len(StripSeq)
                    if SeqLength >= args.min_len:
                        outputfile.write('>%s\n%s\n' % (rec.id, StripSeq))
                    else:
                        errorfile.write('>ERROR:LENGTH=%i|%s\n%s\n' % (SeqLength, orig_id, str(rec.seq)))

def makeDB(input):
    db_details = args.out + '.udb.txt'
    usearch_db = args.out + '.udb'
    Total = amptklib.countfasta(input)
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

#Do a version check
usearch = args.usearch
amptklib.versionDependencyChecks(usearch)

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
RevPrimer = revcomp_lib.RevComp(RevPrimer)
fwdLen = len(FwdPrimer)

if not args.trimming:
    amptklib.log.info("Searching for primers, this may take awhile: Fwd: %s  Rev: %s" % (FwdPrimer, RevPrimer))
else:
    amptklib.log.info("Working on file: %s" % args.fasta)

if not args.cpus:
    cpus = multiprocessing.cpu_count()
else:
    cpus = args.cpus

#create temp directory
pid = os.getpid()
folder = 'amptk_tmp_' + str(pid)
if not os.path.exists(folder):
    os.makedirs(folder)

SeqCount = amptklib.countfasta(args.fasta)
amptklib.log.info('{0:,}'.format(SeqCount) + ' records loaded')
#if only 1 cpu just process data
if cpus == 1:
    stripPrimer(args.fasta)
else:
    amptklib.log.info("Using %i cpus to process data" % cpus)

    #now split it into chunks (as many cpus as are queried)
    amptklib.split_fasta(args.fasta, folder, cpus*2)

    #get list of files
    file_list = []
    for file in os.listdir(folder):
        if file.endswith(".fasta"):
            file = os.path.join(folder, file)
            file_list.append(file)

    #finally process reads over number of cpus
    amptklib.runMultiProgress(stripPrimer, file_list, cpus)

#now concatenate outputs together
OutName = args.out + '.extracted.fa'
ErrorName = args.out + '.errors.fa'
with open(OutName, 'w') as outfile:
    with open(ErrorName, 'w') as outfile2:
        for filename in os.listdir(os.path.join(folder)):
            if filename.endswith('.extracted.fa'):
                if filename == OutName:
                    continue
                with open(os.path.join(folder, filename), 'rU') as readfile:
                    shutil.copyfileobj(readfile, outfile)
            if filename.endswith('.errors.fa'):
                if filename == ErrorName:
                    continue
                with open(os.path.join(folder, filename), 'rU') as readfile:
                    shutil.copyfileobj(readfile, outfile2)

if not args.debug:
    #clean up tmp folder
    shutil.rmtree(folder)

#parse stats from the error fasta file
noID = 0
ambig = 0
tooShort = 0
noPrimer = 0
with open(ErrorName, 'rU') as infile:
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

Passed = amptklib.countfasta(OutName)
amptklib.log.info('{0:,}'.format(Passed) + ' records passed (%.2f%%)' % (Passed*100.0/SeqCount))
amptklib.log.info('Errors: {0:,}'.format(noID) + ' no taxonomy info, ' + 
                  '{0:,}'.format(tooShort) + ' too short, ' +
                  '{0:,}'.format(ambig) + ' too many ambiguous bases, ' +
                  '{0:,}'.format(noPrimer) + ' no primers found')
if args.derep_fulllength:
    amptklib.log.info("Now dereplicating sequences (collapsing identical sequences)")
    derep_tmp = args.out + '.derep.extracted.fa'
    os.rename(OutName, derep_tmp)
    dereplicate(derep_tmp, OutName)
    Total = amptklib.countfasta(OutName)
    amptklib.log.info('{0:,}'.format(Total) + ' records passed (%.2f%%)' % (Total*100.0/SeqCount))
    os.remove(derep_tmp)

if args.create_db:
    makeDB(OutName)

print "-------------------------------------------------------"