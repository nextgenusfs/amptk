import sys, logging, csv, os, subprocess
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

ASCII = {'!':'0','"':'1','#':'2','$':'3','%':'4','&':'5',"'":'6','(':'7',')':'8','*':'9','+':'10',',':'11','-':'12','.':'13','/':'14','0':'15','1':'16','2':'17','3':'18','4':'19','5':'20','6':'21','7':'22','8':'23','9':'24',':':'25',';':'26','<':'27','=':'28','>':'29','?':'30','@':'31','A':'32','B':'33','C':'34','D':'35','E':'36','F':'37','G':'38','H':'39','I':'40','J':'41','K':'42','L':'43','M':'44','N':'45','O':'46','P':'47','Q':'48','R':'49','S':'50'}


class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def myround(x, base=10):
    return int(base * round(float(x)/base))

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count
    
def countfastq(input):
    lines = sum(1 for line in open(input))
    count = int(lines) / 4
    return count

def line_count(fname):
    with open(fname) as f:
        i = -1
        for i, l in enumerate(f):
            pass
    return i + 1

def line_count2(fname):
    count = 0
    with open(fname, 'rU') as f:
        for line in f:
            if not '*' in line:
                count += 1
    return count
    
def getreadlength(input):
    with open(input) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                read_length = len(line) - 1 #offset to switch to 1 based counts
            elif i > 2:
                break
    return read_length
    
def get_vsearch_version():
    version = subprocess.Popen(['vsearch', '--version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    version = version.replace('vsearch v', '')
    version = version.split('_')[0]
    print version
    return version

def checkvsearch():
    vers = get_vsearch_version().split('.')
    if int(vers[0]) > 1:
        return True
    else:
        if int(vers[1]) > 9:
            return True
        else:
            if int(vers[2]) > 0:
                return True
            else:
                return False

def get_version():
    version = subprocess.Popen(['ufits', 'version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    return version

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

def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='[%I:%M:%S %p]')
    else:
        stdoutformat = logging.Formatter(colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='[%I:%M:%S %p]')
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
    
def FastMaxEEFilter(input, trunclen, maxee, output):
    with open(output, 'w') as out:
        with open(input, 'rU') as file:
            for title, seq, qual in FastqGeneralIterator(file):
                trunclen = int(trunclen)
                Seq = seq[:trunclen]
                Qual = qual[:trunclen]
                ee = 0
                for bp, Q in enumerate(Qual):
                    q = int(ASCII.get(Q))
                    P = 10**(float(-q)/10)
                    ee += P
                if ee <= float(maxee):
                    out.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))

def MaxEEFilter(input, maxee):
    with open(input, 'rU') as f:
        for rec in SeqIO.parse(f, "fastq"):
            ee = 0
            for bp, Q in enumerate(rec.letter_annotations["phred_quality"]):
                P = 10**(float(-Q)/10)
                ee += P
            if ee <= float(maxee):
                rec.name = ""
                rec.description = ""
                yield rec
                
def dereplicate(input, output):
    seqs = {}
    with open(input, 'rU') as file:
        for title, sequence, qual in FastqGeneralIterator(file):
            if sequence not in seqs:
                if title.endswith(';'):
                    seqs[sequence]=title+'size=1;'
                else:
                    seqs[sequence]=title+';size=1;'
            else:
                count = int(seqs[sequence].split('=')[-1].rstrip(';')) + 1
                formated_string = seqs[sequence].rsplit('=', 1)[0]+'='+str(count)+';'
                seqs[sequence] = formated_string
    with open(output, 'w') as out:
        for sequence in seqs:
            out.write('>'+seqs[sequence]+'\n'+sequence+'\n')

def convertSize(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Y', suffix) 

def faqual2fastq(fasta, qual, fastq):
    global skipCount
    from Bio.SeqIO.QualityIO import PairedFastaQualIterator
    with open(fastq, 'w') as output:
        records = PairedFastaQualIterator(open(fasta), open(qual))
        for rec in records:
            try:
                SeqIO.write(rec, output, 'fastq')
            except ValueError:
                skipCount +1
    return skipCount
    
def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize
    
def fastqreindex(input, output):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    count = 1
    with open(output, 'w') as out:
        with open(input, 'rU') as fastq:
            for title, sequence, qual in FastqGeneralIterator(fastq):
                cols = title.split(';')
                header = 'R_'+str(count)+';'+cols[1]+';'
                count += 1
                out.write("@%s\n%s\n+\n%s\n" % (header, sequence, qual))


primer_db = {'fITS7': 'GTGARTCATCGAATCTTTG', 'ITS4': 'TCCTCCGCTTATTGATATGC', 'ITS1-F': 'CTTGGTCATTTAGAGGAAGTAA', 'ITS2': 'GCTGCGTTCTTCATCGATGC', 'ITS3': 'GCATCGATGAAGAACGCAGC', 'ITS4-B': 'CAGGAGACTTGTACACGGTCCAG', 'ITS1': 'TCCGTAGGTGAACCTGCGG', 'LR0R': 'ACCCGCTGAACTTAAGC', 'LR2R': 'AAGAACTTTGAAAAGAG', 'JH-LS-369rc': 'CTTCCCTTTCAACAATTTCAC', '16S_V3': 'CCTACGGGNGGCWGCAG', '16S_V4': 'GACTACHVGGGTATCTAATCC', 'ITS3_KYO2': 'GATGAAGAACGYAGYRAA'}

def which(name):
    try:
        with open(os.devnull) as devnull:
            diff = ['tbl2asn', 'dustmasker', 'mafft']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
            else:
                subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)


