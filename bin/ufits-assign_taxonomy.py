#!/usr/bin/env python

import sys, os, re, argparse, logging, subprocess, csv, inspect, multiprocessing
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.primer as primer
import lib.revcomp_lib as revcomp_lib
import lib.ufitslib as ufitslib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

parser=argparse.ArgumentParser(prog='ufits-assign_taxonomy.py', usage="%(prog)s [options] -f <FASTA File>",
    description='''assign taxonomy to OTUs''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--fasta', dest='fasta', required=True, help='FASTA input')
parser.add_argument('-o','--out', dest='out', default='ufits-taxonomy', help='Output file (FASTA)')
parser.add_argument('-m','--method', dest='method', default='hybrid',choices=['utax', 'usearch', 'hybrid', 'rdp', 'blast'], help='Taxonomy method')
parser.add_argument('--utax_db', dest='utax_db', default='ITS2.udb', help='UTAX Reference Database')
parser.add_argument('--utax_cutoff', default=0.8, type=restricted_float, help='UTAX confidence value threshold.')
parser.add_argument('--usearch_db', dest='usearch_db', default='USEARCH.udb', help='USEARCH Reference Database')
parser.add_argument('--usearch_cutoff', default=0.7, type=restricted_float, help='USEARCH percent ID threshold.')
parser.add_argument('-r', '--rdp', dest='rdp', default='/Users/jon/scripts/rdp_classifier_2.10.1/dist/classifier.jar', help='Path to RDP Classifier')
parser.add_argument('--rdp_db', dest='rdp_tax', default='fungalits_unite', choices=['16srrna', 'fungallsu', 'fungalits_warcup', 'fungalits_unite'], help='Training set for RDP Classifier')
parser.add_argument('--rdp_cutoff', default=0.8, type=restricted_float, help='RDP confidence value threshold')
parser.add_argument('--local_blast', help='Path to local Blast DB')
parser.add_argument('-u','--usearch', dest="usearch", default='usearch8', help='USEARCH8 EXE')
parser.add_argument('--append_taxonomy', dest="otu_table", nargs='?', help='Append Taxonomy to OTU table')
parser.add_argument('--only_fungi', action='store_true', help='Retain only fungi in OTU table')
args=parser.parse_args()

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

#remove logfile if exists
log_name = args.out + '.taxonomy.log'
if os.path.isfile(log_name):
    os.remove(log_name)

ufitslib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
ufitslib.log.debug(cmd_args)
print "-------------------------------------------------------"

#initialize script, log system info and usearch version
ufitslib.log.info("Operating system: %s" % sys.platform)


if args.method == 'utax':
    if not args.utax_db:
        ufitslib.log.error("No DB specified, exiting")
        os._exit(1)
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        ufitslib.log.error("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    version = usearch_test.split(" v")[1]
    majorV = version.split(".")[0]
    minorV = version.split(".")[1]
    if int(majorV) < 8 or (int(majorV) >= 8 and int(minorV) < 1):
        ufitslib.log.warning("USEARCH version: %s detected you need v8.1.1756 or above" % usearch_test)
        os._exit(1)
    else:
        ufitslib.log.info("USEARCH version: %s" % usearch_test)
    
    #check for correct DB version
    if not args.utax_db:
        ufitslib.log.error("You must specifiy a UTAX database via --utax_db")
        os._exit(1)
    else: #check if exists
        search_db = os.path.join(parentdir, 'DB', args.utax_db)
        if not os.path.isfile(search_db):
            ufitslib.log.error("%s DB was not found, please run `ufits database` command to create formatted DB" % search_db)
        else:
            utax_db = os.path.join(parentdir, 'DB', args.utax_db)
    #Count records
    ufitslib.log.info("Loading FASTA Records")
    total = countfasta(args.fasta)
    ufitslib.log.info('{0:,}'.format(total) + ' OTUs')
    
    #now run through UTAX
    utax_out = args.out + '.utax.txt'
    ufitslib.log.info("Classifying OTUs with UTAX (USEARCH8)")
    cutoff = str(args.utax_cutoff)
    ufitslib.log.debug("%s -utax %s -db %s -utaxout %s -utax_cutoff %s -strand both" % (usearch, args.fasta, utax_db, utax_out, cutoff))
    subprocess.call([usearch, '-utax', args.fasta, '-db', utax_db, '-utaxout', utax_out, '-utax_cutoff', cutoff, '-strand', 'plus'], stdout = FNULL, stderr = FNULL)
    
    #load results into dictionary for appending to OTU table
    ufitslib.log.debug("Loading results into dictionary")
    with open(utax_out, 'rU') as infile:
        reader = csv.reader(infile, delimiter="\t")
        otuDict = {rows[0]:rows[2] for rows in reader}
    
    
elif args.method == 'usearch':
    if not args.usearch_db:
        ufitslib.log.error("No DB specified, exiting")
        os._exit(1)
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        ufitslib.log.error("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    ufitslib.log.info("USEARCH version: %s" % usearch_test)
    #Count records
    ufitslib.log.info("Loading FASTA Records")
    total = countfasta(args.fasta)
    ufitslib.log.info('{0:,}'.format(total) + ' OTUs')
    
    #check for correct DB version
    if not args.usearch_db:
        ufitslib.log.error("You must specifiy a USEARCH database via --usearch_db")
        os._exit(1)
    else: #check if exists
        search_db = os.path.join(parentdir, 'DB', args.usearch_db)
        if not os.path.isfile(search_db):
            ufitslib.log.error("%s DB was not found, please run `ufits database` command to create formatted DB" % search_db)
        else:
            usearch_db = os.path.join(parentdir, 'DB', args.usearch_db)
    #now run through usearch global
    utax_out = args.out + '.usearch.txt'
    ufitslib.log.info("Global alignment OTUs with usearch_global")
    ufitslib.log.debug("%s -usearch_global %s -db %s -id %s -top_hit_only -output_no_hits -userout %s -userfields query+target+id -strand both" % (usearch, args.fasta, usearch_db, str(args.usearch_cutoff), utax_out))
    subprocess.call([usearch, '-usearch_global', args.fasta, '-db', usearch_db, '-userout', utax_out, '-id', str(args.usearch_cutoff), '-strand', 'both', '-output_no_hits', '-top_hit_only', '-userfields', 'query+target+id'], stdout = FNULL, stderr = FNULL)
    
    #load results into dictionary for appending to OTU table
    ufitslib.log.debug("Loading results into dictionary")
    newTable = []
    with open(utax_out, 'rU') as infile:
        reader = csv.reader(infile, delimiter="\t")
        for line in reader:
            if line[1] == "*":
                newTable.append([line[0], 'No hit'])
            else:
                newTable.append([line[0], re.sub("tax=", "", line[1])+' ('+line[2]+')'])
            
    otuDict = {rows[0]:rows[1] for rows in newTable}
    
    ufitslib.log.info("Done classifying OTUs: %s" % utax_out)

elif args.method == 'hybrid':
    usearch = args.usearch
    try:
        usearch_test = subprocess.Popen([usearch, '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        ufitslib.log.error("%s not found in your PATH, exiting." % usearch)
        os._exit(1)
    version = usearch_test.split(" v")[1]
    majorV = version.split(".")[0]
    minorV = version.split(".")[1]
    if int(majorV) < 8 or (int(majorV) >= 8 and int(minorV) < 1):
        ufitslib.log.warning("USEARCH version: %s detected you need v8.1.1756 or above" % usearch_test)
        os._exit(1)
    else:
        ufitslib.log.info("USEARCH version: %s" % usearch_test)
    
    #check for correct DB version
    if not args.utax_db:
        ufitslib.log.error("You must specifiy a UTAX database via --utax_db")
        os._exit(1)
    else: #check if exists
        search_db = os.path.join(parentdir, 'DB', args.utax_db)
        if not os.path.isfile(search_db):
            ufitslib.log.error("%s DB was not found, please run `ufits database` command to create formatted DB" % search_db)
        else:
            utax_db = os.path.join(parentdir, 'DB', args.utax_db)
    if not args.usearch_db:
        ufitslib.log.error("You must specifiy a USEARCH database via --usearch_db")
        os._exit(1)
    else: #check if exists
        search_db = os.path.join(parentdir, 'DB', args.usearch_db)
        if not os.path.isfile(search_db):
            ufitslib.log.error("%s DB was not found, please run `ufits database` command to create formatted DB" % search_db)
        else:
            usearch_db = os.path.join(parentdir, 'DB', args.usearch_db)
              
    #Count records
    ufitslib.log.info("Loading FASTA Records")
    total = countfasta(args.fasta)
    ufitslib.log.info('{0:,}'.format(total) + ' OTUs')
    
    #now run through UTAX
    utax_out = args.out + '.utax.txt'
    usearch_out = args.out + '.usearch.txt'
    ufitslib.log.info("Classifying OTUs with UTAX (USEARCH8)")
    cutoff = str(args.utax_cutoff)
    ufitslib.log.debug("%s -utax %s -db %s -utaxout %s -utax_cutoff %s -strand both" % (usearch, args.fasta, utax_db, utax_out, cutoff))
    subprocess.call([usearch, '-utax', args.fasta, '-db', utax_db, '-utaxout', utax_out, '-utax_cutoff', cutoff, '-strand', 'plus'], stdout = FNULL, stderr = FNULL)
    
    #load results into dictionary for appending to OTU table
    ufitslib.log.debug("Loading results into dictionary")
    with open(utax_out, 'rU') as infile:
        reader = csv.reader(infile, delimiter="\t")
        utaxDict = {rows[0]:rows[2] for rows in reader}
    
    #run through USEARCH
    ufitslib.log.info("Global alignment OTUs with usearch_global")
    ufitslib.log.debug("%s -usearch_global %s -db %s -id %s -top_hit_only -output_no_hits -userout %s -userfields query+target+id -strand both" % (usearch, args.fasta, usearch_db, str(args.usearch_cutoff), utax_out))
    subprocess.call([usearch, '-usearch_global', args.fasta, '-db', usearch_db, '-userout', usearch_out, '-id', str(args.usearch_cutoff), '-strand', 'both', '-output_no_hits', '-top_hit_only', '-userfields', 'query+target+id'], stdout = FNULL, stderr = FNULL)
    
    #load results into dictionary for appending to OTU table
    ufitslib.log.debug("Loading usearch_global results into dictionary")
    newTable = []
    with open(usearch_out, 'rU') as infile:
        reader = csv.reader(infile, delimiter="\t")
        for line in reader:
            if line[1] == "*":
                utaxLook = utaxDict.get(line[0]) or "No Hit"
                if utaxLook != "No Hit":
                    utaxLook = 'UTAX;' + utaxLook
                newTable.append([line[0], utaxLook])
            elif float(line[2]) < 97:
                utaxLook = utaxDict.get(line[0]) or "No Hit"
                if utaxLook != "No Hit":
                    utaxLook = 'UTAX;' + utaxLook
                newTable.append([line[0], utaxLook])
            else:
                #compare the levels of taxonomy, by counting tax levels
                utaxLook = utaxDict.get(line[0])
                utaxCount = utaxLook.count(',')
                usearchResult = line[1]
                usearchCount = usearchResult.count(',')
                if utaxCount > usearchCount:
                    utaxLook = 'UTAX;' + utaxLook
                    newTable.append([line[0], utaxLook])
                else:
                    newTable.append([line[0], re.sub("tax=", "", usearchResult)])
            
    otuDict = {rows[0]:rows[1] for rows in newTable}

elif args.method == 'rdp':
    #check that classifier is installed
    try:
        rdp_test = subprocess.Popen(['java', '-Xmx2000m', '-jar', args.rdp, 'classify'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        ufitslib.log.error("%s not found in your PATH, exiting." % args.rdp)
        os._exit(1)
    
    #Count records
    ufitslib.log.info("Loading FASTA Records")
    total = countfasta(args.fasta)
    ufitslib.log.info('{0:,}'.format(total) + ' OTUs')
    
    #RDP database
    ufitslib.log.info("Using RDP classifier %s training set" % args.rdp_tax)
    
    #run RDP
    rdp_out = args.out + '.rdp.txt'
    subprocess.call(['java', '-Xmx2000m', '-jar', args.rdp, 'classify', '-g', args.rdp_tax, '-o', rdp_out, '-f', 'fixrank', args.fasta])
    
    #load in results and put into dictionary
    new = []
    removal = ["unidentified", "Incertae", "uncultured", "incertae"]
    remove_exp = [re.compile(x) for x in removal]
    f = csv.reader(open(rdp_out), delimiter='\t')
    for col in f:
        if float(col[19]) > args.rdp_cutoff:
            tax = "k:"+col[2]+" ("+col[4]+"),p:"+col[5]+" ("+col[7]+"),c:"+col[8]+" ("+col[10]+"),o:"+col[11]+" ("+col[13]+"),f:"+col[14]+" ("+col[16]+"),g:"+col[17]+" ("+col[19]+")"
        elif float(col[16]) > args.rdp_cutoff:
            tax = "k:"+col[2]+" ("+col[4]+"),p:"+col[5]+" ("+col[7]+"),c:"+col[8]+" ("+col[10]+"),o:"+col[11]+" ("+col[13]+"),f:"+col[14]+" ("+col[16]+")"
        elif float(col[13]) > args.rdp_cutoff:
            tax = "k:"+col[2]+" ("+col[4]+"),p:"+col[5]+" ("+col[7]+"),c:"+col[8]+" ("+col[10]+"),o:"+col[11]+" ("+col[13]+")"
        elif float(col[10]) > args.rdp_cutoff:
            tax = "k:"+col[2]+" ("+col[4]+"),p:"+col[5]+" ("+col[7]+"),c:"+col[8]+" ("+col[10]+")"
        elif float(col[7]) > args.rdp_cutoff:
            tax = "k:"+col[2]+" ("+col[4]+"),p:"+col[5]+" ("+col[7]+")"
        elif float(col[4]) > args.rdp_cutoff:
            tax = "k:"+col[2]+" ("+col[4]+")"
        else:
            tax = "k:unclassified"
        tax_split = tax.split(",")
        tax = [s for s in tax_split if not any(re.search(s) for re in remove_exp)]
        tax = ",".join(tax)
        line = [col[0],tax]
        new.append(line)
    otuDict = dict(new)

elif args.method == 'blast':
    #check if command line blast installed
    try:
        blast_test = subprocess.Popen(['blastn', '-version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    except OSError:
        ufitslib.log.error("BLASTN not found in your PATH, exiting.")
        os._exit(1)
   
    #Count records
    ufitslib.log.info("Loading FASTA Records")
    total = countfasta(args.fasta)
    ufitslib.log.info('{0:,}'.format(total) + ' OTUs')
    
    #now run blast remotely using NCBI nt database
    blast_out = args.out + '.blast.txt'
    outformat = "6 qseqid sseqid pident stitle"
    if args.local_blast:
        #get number of cpus
        cpus = multiprocessing.cpu_count() - 2
        ufitslib.log.info("Running local BLAST using db: %s" % args.local_blast)
        subprocess.call(['blastn', '-num_threads', str(cpus), '-query', args.fasta, '-db', args.local_blast, '-max_target_seqs', '1', '-outfmt', outformat, '-out', blast_out], stderr = FNULL)
    else:
        ufitslib.log.info("Running BLASTN using NCBI remote nt database, this may take awhile")
        subprocess.call(['blastn', '-query', args.fasta, '-db', 'nt', '-remote', '-max_target_seqs', '1', '-outfmt', outformat, '-out', blast_out], stderr = FNULL)
    
    #load results and reformat
    new = []
    f = csv.reader(open(blast_out), delimiter='\t')
    for col in f:
        query = col[0]
        gbID = col[1].split("|")[3]
        pident = col[2]
        name = col[3]
        tax = gbID + ";" + name + " (" + pident + ")"
        line = [query,tax]
        new.append(line)
    otuDict = dict(new)
    
if not args.otu_table:
    base = args.fasta.split('otus.fa')[0]
    otu_table = base + 'otu_table.txt'
    if not os.path.exists(otu_table):
        otu_table = base + 'otu_table.csv'
        if not os.path.isfile(otu_table):
            otu_table = ""
        else:
            ufitslib.log.info("Found OTU table: %s" % otu_table)
    else:
        ufitslib.log.info("Found OTU table: %s" % otu_table)

else:
    otu_table = args.otu_table

#check if otu_table variable is empty, then load in otu table
if otu_table:
    ufitslib.log.info("Appending taxonomy to OTU table and OTUs")
    end = otu_table.rsplit(".", 1)[1]
    if end == 'txt':
        d = '\t'
    if end == 'csv':
        d = ','
    if not args.out:
        basename = otu_table.rsplit(".", 1)[0]
        taxTable = basename + '.taxonomy.txt'
        otuTax = base + '.otus.taxonomy.fa'
    else:
        taxTable = args.out + '.otu_table.taxonomy.txt'
        otuTax = args.out + '.otus.taxonomy.fa'
    
    #append to OTU table
    counts = 0
    with open(taxTable, 'w') as outTable:
        with open(otu_table, 'rU') as inTable:
            reader = csv.reader(inTable, delimiter=d)
            for line in reader:
                if 'OTUId' in line[0]:
                    line.append('Taxonomy')
                else:
                    tax = otuDict.get(line[0]) or "No Hit"
                    line.append(tax)
                if args.only_fungi:
                    if 'OTUId' in line[0]:
                        join_line = ('\t'.join(str(x) for x in line))
                    else:
                        if 'Fungi' in line[-1]:
                            join_line = ('\t'.join(str(x) for x in line))
                            counts += 1
                        else:
                            continue
                else:
                    join_line = ('\t'.join(str(x) for x in line))
                    counts += 1
                outTable.write("%s\n" % join_line)
    if args.only_fungi:
        nonfungal = total - counts
        ufitslib.log.info("Found %i non-fungal OTUs, writing %i fungal hits to taxonomy OTU table" % (nonfungal, counts))
    #append to OTUs         
    with open(otuTax, 'w') as output:
        with open(args.fasta, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'fasta')
            for rec in SeqRecords:
                tax = otuDict.get(rec.id) or "No hit"
                rec.description = tax
                SeqIO.write(rec, output, 'fasta')
            
    ufitslib.log.info("Taxonomy finished: %s" % taxTable)
else:
    ufitslib.log.info("Unable to automatically detect OTU table, skipping append taxonomy.  Try to specifiy path to OTU table.")

print "-------------------------------------------------------"
