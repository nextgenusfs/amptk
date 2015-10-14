#!/usr/bin/env python

import sys, os, re, argparse, logging, subprocess, inspect
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.progress as progress

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-extract_region.py', usage="%(prog)s [options] -f <FASTA File>",
    description='''Script searches for primers and removes them if found.  Useful for trimming a reference dataset for assigning taxonomy after OTU clustering.  It is also capable of reformatting UNITE taxonomy fasta headers to be compatible with UTAX and creating USEARCH/UTAX UBD databases for assigning taxonomy.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fasta', dest='fasta', required=True, help='FASTA input')
parser.add_argument('-o','--out', dest='out', default='ufits', help='Base Name Output files')
parser.add_argument('-f','--fwd_primer', dest='fwd_primer', default='GTGARTCATCGAATCTTTG', help='Forward primer (fITS7)')
parser.add_argument('-r','--rev_primer', dest='rev_primer', default='TCCTCCGCTTATTGATATGC', help='Reverse primer (ITS4)')
parser.add_argument('--skip_trimming', dest='trimming', action='store_true', help='Skip Primer trimming (not recommended)')
parser.add_argument('--unite2utax', dest='utax', default='on', choices=['on', 'off'], help='Reformat UNITE FASTA headers for UTAX')
parser.add_argument('--drop_ns', dest='drop_ns', type=int, default=8, help="Drop Seqeunces with more than X # of N's")
parser.add_argument('--create_db', dest='create_db', choices=['utax', 'usearch'], help="Create USEARCH DB")
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
args=parser.parse_args()

def stripPrimer(records):
    global SeqCount, OutCount, RevCount, NoMatch	  
    for rec in records:
        if SeqCount == 0:
		    progress.InitFile(input)
    
        progress.File2("%u seqs processed, %u seqs pass, %u rev_comp'd pass, %u seqs skipped" % \
        (SeqCount, OutCount, RevCount, NoMatch))

        SeqCount += 1
        if args.utax == 'on':
            fields = rec.description.split("|")
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
            p = tf[1]
            c = tf[2]
            o = tf[3]
            f = tf[4]
            g = tf[5]
            s = tf[6]
            reformat_tax = []
            removal = ("unidentified", "Incertae", "uncultured")
            sp_removal = (" sp", "_sp", "uncultured")
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
            if rec.id.endswith(";tax="): #if there is no taxonomy, get rid of it
                rec.id = ""
            rec.name = ""
            rec.description = ""
        if not args.trimming:
            Seq = rec.seq
            MAX_PRIMER_MISMATCHES = 2
            RevPrimer = revcomp_lib.RevComp(rev_primer)
            BestPosFor, BestDiffsFor = primer.BestMatch2(Seq, fwd_primer, MAX_PRIMER_MISMATCHES)
            if BestDiffsFor < MAX_PRIMER_MISMATCHES:
                if BestPosFor > 0:
                    stripfwdlen = fwdLen + BestPosFor
                    StripSeq = Seq[stripfwdlen:]
                
                    #now look for reverse
                    BestPosRev, BestDiffsRev = primer.BestMatch2(StripSeq, RevPrimer, MAX_PRIMER_MISMATCHES)
                    if BestDiffsRev < MAX_PRIMER_MISMATCHES:
                        StrippedSeq = StripSeq[:BestPosRev]
                    else:
                        StrippedSeq = StripSeq
                    #after stripping primers, check for ambig bases
                    if args.drop_ns != 0 and 'N'*args.drop_ns in StrippedSeq:
                        continue
                    rec.seq = StrippedSeq
                    if rec.seq and rec.id:
                        OutCount += 1
                        yield rec
            else: #if can't find forward primer, try to reverse complement and look again
                RevSeq = revcomp_lib.RevComp(Seq)
                BestPosFor, BestDiffsFor = primer.BestMatch2(RevSeq, fwd_primer, MAX_PRIMER_MISMATCHES)
                if BestDiffsFor < MAX_PRIMER_MISMATCHES:
                    if BestPosFor > 0:
                        stripfwdlen = fwdLen + BestPosFor
                        StripSeq = Seq[stripfwdlen:]
                
                        #now look for reverse
                        BestPosRev, BestDiffsRev = primer.BestMatch2(StripSeq, RevPrimer, MAX_PRIMER_MISMATCHES)
                        if BestDiffsRev < MAX_PRIMER_MISMATCHES:
                            StrippedSeq = StripSeq[:BestPosRev]
                        else:
                            StrippedSeq = StripSeq
                        #after stripping primers, check for ambig bases
                        if args.drop_ns != 0 and 'N'*args.drop_ns in StrippedSeq:
                            continue
                        rec.seq = StrippedSeq
                        if rec.seq and rec.id:
                            RevCount += 1
                            yield rec
                else:
                    NoMatch += 1
        else:
            yield rec

def makeDB(input):
    #need usearch for this, test to make sure version is ok with utax
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        log.error("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    version = usearch_test.split(" v")[1]
    majorV = version.split(".")[0]
    minorV = version.split(".")[1]
    if int(majorV) < 8 or (int(majorV) >= 8 and int(minorV) < 1):
        log.warning("USEARCH version: %s detected you need v8.1.1756 or above" % usearch_test)
        os._exit(1)
    else:
        log.info("USEARCH version: %s" % usearch_test)
    
    db_details = args.out + '.' + args.create_db + '.udb.txt'
    usearch_db = args.out + '.' + args.create_db +'.udb'
    if args.trimming:
        args.fwd_primer = 'None'
        args.rev_primer = 'None'
    db_string = args.fasta + ' ' + args.fwd_primer + ' ' + args.rev_primer + ' ' + str(Total)
    with open(db_details, 'w') as details:
        details.write(db_string)
    report = args.out + '.report.txt'
    
    
    if args.create_db == 'utax':
        #create log file for this to troubleshoot
        utax_log = args.out + '.utax.log'
        if os.path.isfile(utax_log):
            os.remove(utax_log)    
        utaxLog = open(utax_log, 'w')  
        log.info("Creating UTAX Database, this may take awhile")
        log.debug("%s -makeudb_utax %s -output %s -report %s -utax_trainlevels kpcofgs -utax_splitlevels NVpcofgs" % (usearch, input, usearch_db, report))
        subprocess.call([usearch, '-makeudb_utax', input, '-output', usearch_db, '-report', report, '-utax_trainlevels', 'kpcofgs', '-utax_splitlevels', 'NVpcofgs', '-notrunclabels'], stdout = utaxLog, stderr = utaxLog)
        utaxLog.close()
        log.info("Database %s created successfully" % usearch_db)
        
    if args.create_db == 'usearch':
        log.info("Creating USEARCH Database")
        log.debug("%s -makeudb_usearch %s -output %s -notrunclabels" % (usearch, input, usearch_db))
        subprocess.call([usearch, '-makeudb_usearch', input, '-output', usearch_db, '-notrunclabels'], stdout = FNULL, stderr = FNULL)
        log.info("Database %s created successfully" % usearch_db)

def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    else:
        stdoutformat = logging.Formatter(col.GRN+'%(asctime)s'+col.END+': %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)

#remove logfile if exists
log_name = args.out + '.log'
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print "-------------------------------------------------------"

SeqCount = 0
OutCount = 0
RevCount = 0
NoMatch = 0

FileName = args.fasta
fwd_primer = args.fwd_primer
rev_primer = args.rev_primer
fwdLen = len(fwd_primer)

if not args.trimming:
    log.info("Searching for primers, this may take awhile")
else:
    log.info("Working on file: %s" % FileName)

OutName = args.out + '.extracted.fa'
with open(OutName, 'w') as output:
    with open(FileName, 'rU') as input:
        SeqRecords = SeqIO.parse(FileName, 'fasta')
        SeqIO.write(stripPrimer(SeqRecords), output, 'fasta')
progress.FileDone("%u seqs processed, %u seqs pass, %u rev_comp'd pass, %u seqs skipped" % \
        (SeqCount, OutCount, RevCount, NoMatch))
Total = OutCount + RevCount
if args.trimming:
    Total = SeqCount    
log.info("Stats for extracting region: \
\n%10u seqs input \
\n%10u fwd primer found \
\n%10u fwd primer found in revcomp'd seq \
\n%10u total seqs output (%.2f%%)" % (SeqCount, OutCount, RevCount, Total, Total*100.0/SeqCount))

if args.create_db:
    makeDB(OutName)

print "-------------------------------------------------------"