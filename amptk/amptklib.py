from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import logging
import csv
import os
import subprocess
import multiprocessing
import platform
import time
import shutil
import gzip
import edlib
import ast
from builtins import range
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from natsort import natsorted
from amptk.__version__ import __version__
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
    
parentdir = os.path.join(os.path.dirname(__file__))

ASCII = {'!':'0','"':'1','#':'2','$':'3','%':'4','&':'5',
         "'":'6','(':'7',')':'8','*':'9','+':'10',',':'11',
         '-':'12','.':'13','/':'14','0':'15','1':'16','2':'17',
         '3':'18','4':'19','5':'20','6':'21','7':'22','8':'23',
         '9':'24',':':'25',';':'26','<':'27','=':'28','>':'29',
         '?':'30','@':'31','A':'32','B':'33','C':'34','D':'35',
         'E':'36','F':'37','G':'38','H':'39','I':'40','J':'41',
         'K':'42','L':'43','M':'44','N':'45','O':'46','P':'47',
         'Q':'48','R':'49','S':'50'}

primer_db = {'fITS7': 'GTGARTCATCGAATCTTTG',
            'fITS7-ion': 'AGTGARTCATCGAATCTTTG',
            'ITS4': 'TCCTCCGCTTATTGATATGC', 
            'ITS1-F': 'CTTGGTCATTTAGAGGAAGTAA', 
            'ITS2': 'GCTGCGTTCTTCATCGATGC', 
            'ITS3': 'GCATCGATGAAGAACGCAGC', 
            'ITS4-B': 'CAGGAGACTTGTACACGGTCCAG', 
            'ITS1': 'TCCGTAGGTGAACCTGCGG', 
            'LR0R': 'ACCCGCTGAACTTAAGC', 
            'LR2R': 'AAGAACTTTGAAAAGAG', 
            'JH-LS-369rc': 'CTTCCCTTTCAACAATTTCAC', 
            '16S_V3': 'CCTACGGGNGGCWGCAG', 
            '16S_V4': 'GACTACHVGGGTATCTAATCC', 
            'ITS3_KYO2': 'GATGAAGAACGYAGYRAA', 
            'COI-F': 'GGTCAACAAATCATAAAGATATTGG', 
            'COI-R': 'GGWACTAATCAATTTCCAAATCC', 
            '515FB': 'GTGYCAGCMGCCGCGGTAA', 
            '806RB': 'GGACTACNVGGGTWTCTAAT', 
            'ITS4-B21': 'CAGGAGACTTGTACACGGTCC',
            'LCO1490': 'GGTCAACAAATCATAAAGATATTGG',
            'mlCOIintR': 'GGRGGRTASACSGTTCASCCSGTSCC'}


degenNuc = [("R", "A"), ("R", "G"), 
            ("M", "A"), ("M", "C"),
            ("W", "A"), ("W", "T"),
            ("S", "C"), ("S", "G"),
            ("Y", "C"), ("Y", "T"),
            ("K", "G"), ("K", "T"),
            ("V", "A"), ("V", "C"), ("V", "G"),
            ("H", "A"), ("H", "C"), ("H", "T"),
            ("D", "A"), ("D", "G"), ("D", "T"),
            ("B", "C"), ("B", "G"), ("B", "T"),
            ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C"),
            ("X", "G"), ("X", "A"), ("X", "T"), ("X", "C")]

degenNucSimple = [("R", "A"), ("R", "G"), 
            ("M", "A"), ("M", "C"),
            ("W", "A"), ("W", "T"),
            ("S", "C"), ("S", "G"),
            ("Y", "C"), ("Y", "T"),
            ("K", "G"), ("K", "T"),
            ("V", "A"), ("V", "C"), ("V", "G"),
            ("H", "A"), ("H", "C"), ("H", "T"),
            ("D", "A"), ("D", "G"), ("D", "T"),
            ("B", "C"), ("B", "G"), ("B", "T")]

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

#setup functions
def download(url, name):
    file_name = name
    u = urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = 0
    for x in meta.items():
        if x[0].lower() == 'content-length':
            file_size = int(x[1])
    print("Downloading: {0} Bytes: {1}".format(url, file_size))
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        file_size_dl += len(buffer)
        f.write(buffer)
        p = float(file_size_dl) / file_size
        status = r"{0}  [{1:.2%}]".format(file_size_dl, p)
        status = status + chr(8)*(len(status)+1)
        sys.stdout.write(status)
    f.close()

#functions for system checks, etc
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = open(os.devnull, 'w'), cwd=parentdir).communicate()[0]
        return out
    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', '--short', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = False
    return GIT_REVISION

def getSize(filename):
    st = os.stat(filename)
    return st.st_size

def checkfile(input):
    if os.path.isfile(input):
        filesize = getSize(input)
        if int(filesize) < 1:
            return False
        else:
            return True
    elif os.path.islink(input):
        return True
    else:
        return False

def SafeRemove(input):
    if os.path.isdir(input):
        shutil.rmtree(input)
    elif os.path.isfile(input):
        os.remove(input)
    else:
        return

def number_present(s):
    return any(i.isdigit() for i in s)
    
def get_version():
    if git_version():
        version = __version__+'-'+git_version
    else:
        version = __version__
    return version

def get_usearch_version(usearch):
    try:
        version = subprocess.Popen([usearch, '-version'], stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    except OSError:
        log.error("%s not found in your PATH, exiting." % usearch)
        sys.exit(1)
    vers = version[0].decode('utf-8').replace('usearch v', '')
    vers2 = vers.split('_')[0]
    return vers2

def get_vsearch_version():
    version = subprocess.Popen(['vsearch', '--version'], stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    for v in version:
        v = v.decode('utf-8')
        if v.startswith('vsearch'):
            vers = v.replace('vsearch v', '')
            vers2 = vers.split('_')[0]
            return vers2

def versiontuple(v):
    return tuple(map(int, (v.split("."))))

def gvc(input, check):
    if versiontuple(input) >= versiontuple(check):
        return True
    else:
        return False

def versionDependencyChecks(usearch):
    #to run amptk need usearch > 9.0.2132 and vsearch > 2.2.0
    amptk_version = get_version()
    usearch_version = get_usearch_version(usearch)
    vsearch_version = get_vsearch_version()
    usearch_pass = '9.0.2132'
    vsearch_pass = '2.2.0'
    if not gvc(usearch_version, usearch_pass):
        log.error("USEARCH v%s detected, needs to be atleast v%s" % (usearch_version, usearch_pass))
        sys.exit(1)
    if not gvc(vsearch_version, vsearch_pass):
        log.error("VSEARCH v%s detected, needs to be atleast v%s" % (vsearch_version, vsearch_pass))
        sys.exit(1)
    log.info("AMPtk v%s, USEARCH v%s, VSEARCH v%s" % (amptk_version, usearch_version, vsearch_version))

def checkusearch10(usearch):
    #to run amptk need usearch > 10.0.2132 and vsearch > 2.2.0
    amptk_version = get_version()
    usearch_version = get_usearch_version(usearch)
    vsearch_version = get_vsearch_version()
    usearch_pass = '10.0.240'
    vsearch_pass = '2.2.0'
    if not gvc(usearch_version, usearch_pass):
        log.error("USEARCH v%s detected, needs to be atleast v%s" % (usearch_version, usearch_pass))
        sys.exit(1)
    if not gvc(vsearch_version, vsearch_pass):
        log.error("VSEARCH v%s detected, needs to be atleast v%s" % (vsearch_version, vsearch_pass))
        sys.exit(1)
    log.info("AMPtk v%s, USEARCH v%s, VSEARCH v%s" % (amptk_version, usearch_version, vsearch_version))

def checkRversion():
    #need to have R version > 3.2
    cmd = ['Rscript', '--vanilla', os.path.join(parentdir, 'check_version.R')]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    versions = stdout.decode('utf-8').replace(' \n', '')
    Rvers = versions.split(',')[0]
    dada2 = versions.split(',')[1]
    phyloseq = versions.split(',')[2]
    lulu = versions.split(',')[3]
    return (Rvers, dada2, phyloseq, lulu)

def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize

def getCPUS():
    cores = multiprocessing.cpu_count()
    return cores

def MemoryCheck():
    import psutil
    mem = psutil.virtual_memory()
    RAM = int(mem.total)
    return round(RAM / 1024000000)
                   
def systemOS():
    if sys.platform == 'darwin':
        system_os = 'MacOSX '+ platform.mac_ver()[0]
    elif sys.platform == 'linux':
        linux_version = platform.linux_distribution()
        system_os = linux_version[0]+ ' '+linux_version[1]
    else:
        system_os = sys.platform
    return system_os

def SystemInfo():
    system_os = systemOS()
    python_vers = str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])   
    log.info("OS: %s, %i cores, ~ %i GB RAM. Python: %s" % (system_os, multiprocessing.cpu_count(), MemoryCheck(), python_vers))
    #print python modules to logfile
    mod_versions()
    
def mod_versions():
    import pkg_resources
    modules = ['numpy', 'pandas', 'matplotlib', 'psutil', 'natsort', 'biopython', 'edlib', 'biom-format']
    results = []
    for x in modules:
        try:
            vers = pkg_resources.get_distribution(x).version
            hit = "%s v%s" % (x, vers)
        except pkg_resources.DistributionNotFound:
            hit = "%s NOT installed!" % x
        results.append(hit)
        if x == 'edlib':
            if '.post' in vers:
                vers = vers.split('.post')[0]
            if not gvc(vers, '1.2.1'):
                log.error("Edlib v%s detected, at least v1.2.1 required for degenerate nucleotide search, please upgrade.  e.g. pip install -U edlib or conda install edlib" % vers)
                sys.exit(1)
    log.debug("Python Modules: %s" % ', '.join(results))
 
  
class gzopen(object):
    def __init__(self, fname):
        f = open(fname)
        # Read magic number (the first 2 bytes) and rewind.
        magic_number = f.read(2)
        f.seek(0)
        # Encapsulated 'self.f' is a file or a GzipFile.
        if magic_number == '\x1f\x8b':
            self.f = gzip.GzipFile(fileobj=f)
        else:
            self.f = f
    # Define '__enter__' and '__exit__' to use in
    # 'with' blocks. Always close the file and the
    # GzipFile if applicable.
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        try:
            self.f.fileobj.close()
        except AttributeError:
            pass
        finally:
            self.f.close()
    # Reproduce the interface of an open file
    # by encapsulation.
    def __getattr__(self, name):
        return getattr(self.f, name)
    def __iter__(self):
        return iter(self.f)
    def __next__(self):
        return next(self.f)
      
def Funzip(input, output, cpus):
    '''
    function to unzip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '--decompress', '-c', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '--decompress', '-c', input]
    try:
        runSubprocess2(cmd, log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)

def Fzip(input, output, cpus):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '-c', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '-c', input]
    try:
        runSubprocess2(cmd, log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)

def Fzip_inplace(input):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    cpus = multiprocessing.cpu_count()
    if which('pigz'):
        cmd = ['pigz', '-f', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '-f', input]
    try:
        runSubprocess(cmd, log)
    except NameError:
        subprocess.call(cmd)

def fastalen2dict(input):
    Lengths = {}
    with open(input, 'r') as infile:
        for rec in SeqIO.parse(infile, 'fasta'):
            if not rec.id in Lengths:
                Lengths[rec.id] = len(rec.seq)
    return Lengths

def myround(x, base=10):
    return int(base * round(float(x)/base))

def GuessRL(input):
    #read first 50 records, get length then exit
    lengths = []
    for title, seq, qual in FastqGeneralIterator(gzopen(input)):
        if len(lengths) < 50:
            lengths.append(len(seq))
        else:
            break
    return myround(max(set(lengths)))

def countfasta(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count
    
def countfastq(input):
    lines = sum(1 for line in gzopen(input))
    count = int(lines) // 4
    return count

def line_count(fname):
    with open(fname) as f:
        i = -1
        for i, l in enumerate(f):
            pass
    return i + 1
    
def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

def line_count2(fname):
    count = 0
    with open(fname, 'r') as f:
        for line in f:
            if not '*' in line:
                count += 1
    return count
    
def getreadlength(input):
    with gzopen(input) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                read_length = len(line) - 1 #offset to switch to 1 based counts
            elif i > 2:
                break
    return read_length
    
def runSubprocess(cmd, logfile):
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout:
        logfile.debug(stdout)
    if stderr:
        logfile.debug(stderr)

def runSubprocess2(cmd, logfile, output):
    #function where output of cmd is STDOUT, capture STDERR in logfile
    logfile.debug(' '.join(cmd))
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] != None:
            logfile.debug(stderr)

def runSubprocess3(cmd, logfile, folder, output):
    #function where output of cmd is STDOUT, capture STDERR in logfile
    logfile.debug(' '.join(cmd))
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=out, cwd=folder)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] != None:
            logfile.debug(stderr)
            
def runSubprocess4(cmd, logfile, logfile2):
    #function where cmd is issued in logfile, and log captured in logfile 2
    logfile.debug(' '.join(cmd))
    with open(logfile2, 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=out)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] != None:
            logfile.debug(stderr)

def runSubprocess5(cmd):
    #function where no logfile and stdout/stderr to fnull
    FNULL = open(os.devnull, 'w')
    #print(' '.join(cmd))
    subprocess.call(cmd, stdout=FNULL, stderr=FNULL)

def getSize(filename):
    st = os.stat(filename)
    return st.st_size

def check_valid_file(input):
    if os.path.isfile(input):
        filesize = getSize(input)
        if int(filesize) < 1:
            return False
        else:
            return True
    else:
        return False
        
def bam2fastq(input, output):
    from . import pybam
    with open(output, 'w') as fastqout:
        with open(input, 'r') as bamin:
            for title, seq, qual in pybam.read(bamin,['sam_qname', 'sam_seq','sam_qual']):
                fastqout.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with gzopen(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def return_lines(path, linepos, nstart, nstop):
    """return nsamp lines from path where line offsets are in linepos"""
    offsets = linepos[int(nstart):int(nstop)]
    lines = []
    with gzopen(path) as inf:
        for offset in offsets:
            inf.seek(offset)
            lines.append(inf.readline())
    return lines     

def split_fastq(input, numseqs, outputdir, chunks):
    #get number of sequences and then number of sequences in each chunk
    numlines = numseqs*4
    n = numlines // chunks
    #make sure n is divisible by 4 (number of lines in fastq)
    if ( n % 4 ) != 0:
        n = ((n // 4) + 1) * 4
    splits = []
    count = 0
    for i in range(chunks):
        start = count
        end = count+n
        if end > numlines:
            end = numlines
        splits.append((start, end))
        count += n
    #get line positions from file
    linepos = scan_linepos(input)
    #make sure output directory exists
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    #loop through the positions and write output
    for i, x in enumerate(splits):
        num = i+1
        with open(os.path.join(outputdir, 'chunk_'+str(num)+'.fq'), 'w') as output:
            lines = return_lines(input, linepos, x[0], x[1])
            output.write('%s' % ''.join(lines))
            
def split_fastqPE(R1, R2, numseqs, outputdir, chunks):
    #get number of sequences and then number of sequences in each chunk
    numlines = numseqs*4
    n = numlines // chunks
    #make sure n is divisible by 4 (number of lines in fastq)
    if ( n % 4 ) != 0:
        n = ((n // 4) + 1) * 4
    splits = []
    count = 0
    for i in range(chunks):
        start = count
        end = count+n
        if end > numlines:
            end = numlines
        splits.append((start, end))
        count += n
    #get line positions from file
    linepos = scan_linepos(R1)
    #make sure output directory exists
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    #loop through the positions and write output
    for i, x in enumerate(splits):
        num = i+1
        with open(os.path.join(outputdir, 'chunk_'+str(num)+'_R1.fq'), 'w') as output1:
            with open(os.path.join(outputdir, 'chunk_'+str(num)+'_R2.fq'), 'w') as output2:
                lines1 = return_lines(R1, linepos, x[0], x[1])
                output1.write('%s' % ''.join(lines1))
                lines2 = return_lines(R2, linepos, x[0], x[1])
                output2.write('%s' % ''.join(lines2))   

def split_fastqPEandI(R1, R2, I1, numseqs, outputdir, chunks):
    #get number of sequences and then number of sequences in each chunk
    numlines = numseqs*4
    n = numlines // chunks
    #make sure n is divisible by 4 (number of lines in fastq)
    if ( n % 4 ) != 0:
        n = ((n // 4) + 1) * 4
    splits = []
    count = 0
    for i in range(chunks):
        start = count
        end = count+n
        if end > numlines:
            end = numlines
        splits.append((start, end))
        count += n
    #get line positions from file
    linepos = scan_linepos(R1)
    linepos2 = scan_linepos(I1)
    #make sure output directory exists
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    #loop through the positions and write output
    for i, x in enumerate(splits):
        num = i+1
        with open(os.path.join(outputdir, 'chunk_'+str(num)+'_R1.fq'), 'w') as output1:
            with open(os.path.join(outputdir, 'chunk_'+str(num)+'_R2.fq'), 'w') as output2:
                with open(os.path.join(outputdir, 'chunk_'+str(num)+'_R3.fq'), 'w') as output3:
                    lines1 = return_lines(R1, linepos, x[0], x[1])
                    output1.write('%s' % ''.join(lines1))
                    lines2 = return_lines(R2, linepos, x[0], x[1])
                    output2.write('%s' % ''.join(lines2))
                    lines3 = return_lines(I1, linepos2, x[0], x[1])
                    output3.write('%s' % ''.join(lines3))

def split_fasta(input, outputdir, chunks):
    #function to return line positions of fasta files for chunking
    fastapos = []
    position = 0
    numseqs = 0
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                numseqs += 1
                fastapos.append(position)
            position += 1
    splits = []
    n = numseqs // chunks
    num = 0
    lastline = line_count(input)
    for i in range(chunks):
        if i == 0:
            start = 0
            num = n
            lastpos = fastapos[n+1]
        elif i == chunks-1: #last one
            start = lastpos
            lastpos = lastline
        else:
            start = lastpos
            num = num + n
            try:
                lastpos = fastapos[num+1] #find the n+1 seq
            except IndexError:
                lastpos = fastapos[-1]
        splits.append((start, lastpos))
    #check if output folder exists, if not create it
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    #get line positions from file
    linepos = scan_linepos(input)
    #loop through the positions and write output
    for i, x in enumerate(splits):
        num = i+1
        with open(os.path.join(outputdir, 'chunk_'+str(num)+'.fasta'), 'w') as output:
            lines = return_lines(input, linepos, x[0], x[1])
            output.write('{:}'.format(''.join(lines)))

def trim3prime(input, trimlen, output, removelist):
    with open(output, 'w') as outfile:
        for title, seq, qual in FastqGeneralIterator(gzopen(input)):
            if not title.split(' ')[0] in removelist:
                Seq = seq[:trimlen]
                Qual = qual[:trimlen]
                outfile.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))

def PEsanitycheck(R1, R2):
    R1count = line_count(R1)
    R2count = line_count(R2)
    if R1count != R2count:
        return False
    else:
        return True

def PEandIndexCheck(R1, R2, R3):
    R1count = countfastq(R1)
    R2count = countfastq(R2)
    R3count = countfastq(R3)
    if R1count == R2count == R3count:
        return True
    else:
        return False

def mapIndex(seq, mapDict, bcmismatch):
    besthit = []
    for index_name,index in mapDict.items():
        align = edlib.align(index, seq, mode="SHW", k=bcmismatch, additionalEqualities=degenNuc)
        if align["editDistance"] < 0:
            continue
        elif align["editDistance"] == 0:
            return (index_name, 0)
        else:
            if len(besthit) < 3:
                besthit = [index, index_name, align["editDistance"]]
            else:
                if align["editDistance"] < int(besthit[2]):
                    besthit = [index, index_name, align["editDistance"]]
    if len(besthit) == 3:
        return (besthit[1], besthit[2])
    else:
        return (None,None)
 
    
def DemuxIllumina(R1, R2, I1, mapDict, mismatch, fwdprimer, revprimer, primer_mismatch, outR1, outR2):
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest
    Total = 0
    FPrimer = 0
    RPrimer = 0
    BCFound = 0
    #function to loop through PE reads, renaming according to index
    file1 = FastqGeneralIterator(gzopen(R1))
    file2 = FastqGeneralIterator(gzopen(R2))
    file3 = FastqGeneralIterator(gzopen(I1))
    counter = 1
    with open(outR1, 'w') as outfile1:
        with open(outR2, 'w') as outfile2:
            for read1, read2, index in zip(file1, file2, file3):
                Total += 1
                Name,Diffs = mapIndex(index[1], mapDict, mismatch)
                if Name:
                    BCFound += 1
                    #strip primers if found
                    R1ForPos = trimForPrimer(fwdprimer, read1[1], primer_mismatch)
                    R1RevPos = trimRevPrimer(revprimer, read1[1], primer_mismatch)
                    R2ForPos = trimForPrimer(revprimer, read2[1], primer_mismatch)
                    R2RevPos = trimRevPrimer(fwdprimer, read2[1], primer_mismatch)
                    if R1ForPos > 0:
                        FPrimer += 1
                    if R1RevPos > 0:
                        RPrimer += 1
                    header = 'R_'+str(counter)+';barcodelabel='+Name+';bcseq='+index[1]+';bcdiffs='+str(Diffs)+';'
                    outfile1.write('@%s\n%s\n+\n%s\n' % (header, read1[1][R1ForPos:R1RevPos], read1[2][R1ForPos:R1RevPos]))
                    outfile2.write('@%s\n%s\n+\n%s\n' % (header, read2[1][R2ForPos:R2RevPos], read2[2][R2ForPos:R2RevPos]))
                    counter += 1
    return Total, BCFound, FPrimer, RPrimer
                    
def stripPrimersPE(R1, R2, RL, samplename, fwdprimer, revprimer, primer_mismatch, require_primer, full_length, outR1, outR2):
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest
    #can walk through dataset in pairs
    file1 = FastqGeneralIterator(gzopen(R1))
    file2 = FastqGeneralIterator(gzopen(R2))
    counter = 1
    Total = 0
    multihits = 0
    findForPrimer = 0
    findRevPrimer = 0
    with open(outR1, 'w') as outfile1:
        with open(outR2, 'w') as outfile2:
            for read1, read2 in zip(file1, file2):
                Total += 1
                ffp = False
                frp = False
                R1Seq = read1[1][:RL]
                R1Qual = read1[2][:RL]
                R2Seq = read2[1][:RL]
                R2Qual = read2[2][:RL]
                ForTrim, RevTrim = (0,)*2
                #look for forward primer in forward read
                R1foralign = edlib.align(fwdprimer, R1Seq, mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)            
                if R1foralign['editDistance'] < 0:
                    if require_primer == 'on' or full_length: #not found
                        continue
                else:
                    if len(R1foralign['locations']) > 1: #multiple hits
                        multihits += 1
                        continue
                    try:
                        ForTrim = R1foralign["locations"][0][1]+1
                        findForPrimer += 1
                        ffp = True
                    except IndexError:
                        pass
                R1revalign = edlib.align(RevComp(revprimer), R1Seq, mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                if R1revalign['editDistance'] < 0:
                    R1RevCut = RL
                else:
                    R1RevCut = R1revalign["locations"][0][0]
                    findRevPrimer += 1
                    frp = True
                #look for reverse primer in reverse read
                R2foralign = edlib.align(revprimer, R2Seq, mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                if R2foralign['editDistance'] < 0:
                    if require_primer == 'on' or full_length: #not found
                        continue
                else:
                    if len(R2foralign['locations']) > 1: #multiple hits
                        multihits += 1
                        continue
                    try:
                        RevTrim = R2foralign["locations"][0][1]+1
                        if not frp:
                            findRevPrimer += 1
                    except IndexError:
                        pass        
                R2revalign = edlib.align(RevComp(fwdprimer), R2Seq, mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                if R2revalign['editDistance'] < 0:
                    R2RevCut = RL
                else:
                    R2RevCut = R2revalign["locations"][0][0]
                    if not ffp:
                        findForPrimer += 1               
                header = 'R_{:};barcodelabel={:};'.format(counter,samplename)
                outfile1.write('@%s\n%s\n+\n%s\n' % (header, R1Seq[ForTrim:R1RevCut], R1Qual[ForTrim:R1RevCut]))
                outfile2.write('@%s\n%s\n+\n%s\n' % (header, R2Seq[RevTrim:R2RevCut], R2Qual[RevTrim:R2RevCut]))
                counter += 1
    return Total, counter-1, multihits, findForPrimer, findRevPrimer

def primerFound(primer, seq, mismatch):
    align = edlib.align(primer, seq, mode="HW", k=mismatch, additionalEqualities=degenNuc)
    if align['editDistance'] < 0:
        return False
    else:
        return True

def illuminaReorient(R1, R2, fwdprimer, revprimer, mismatch, ReadLength, outR1, outR2):
    '''
    function to re-orient reads based on primer sequences
    only useful if primers in the reads
    drops reads that don't have matching primers
    '''
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest
    Total = 0
    Correct = 0
    Flipped = 0
    Dropped = 0
    #function to loop through PE reads, renaming according to index
    file1 = FastqGeneralIterator(gzopen(R1))
    file2 = FastqGeneralIterator(gzopen(R2))
    with open(outR1, 'w') as outfile1:
        with open(outR2, 'w') as outfile2:
            for read1, read2 in zip(file1, file2):
                Total += 1
                if primerFound(fwdprimer, read1[1], mismatch) and primerFound(revprimer, read2[1], mismatch):
                    Correct += 1
                    outfile1.write('@%s\n%s\n+\n%s\n' % (read1[0], read1[1][:ReadLength], read1[2][:ReadLength]))
                    outfile2.write('@%s\n%s\n+\n%s\n' % (read2[0], read2[1][:ReadLength], read2[2][:ReadLength]))
                elif primerFound(fwdprimer, read2[1], mismatch) and primerFound(revprimer, read1[1], mismatch):
                    Flipped += 1
                    outfile1.write('@%s\n%s\n+\n%s\n' % (read2[0], read2[1][:ReadLength], read2[2][:ReadLength]))
                    outfile2.write('@%s\n%s\n+\n%s\n' % (read1[0], read1[1][:ReadLength], read1[2][:ReadLength]))
                else:
                    Dropped += 1
    return Total, Correct, Flipped, Dropped
        

def demuxIlluminaPE(R1, R2, fwdprimer, revprimer, samples, forbarcodes, revbarcodes, barcode_mismatch, primer_mismatch, outR1, outR2, stats):
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest
    #function to loop through PE reads, renaming according to index
    file1 = FastqGeneralIterator(gzopen(R1))
    file2 = FastqGeneralIterator(gzopen(R2))
    counter = 1
    Total = 0
    NoBarcode = 0
    NoRevBarcode = 0
    NoPrimer = 0
    NoRevPrimer = 0
    ValidSeqs = 0
    with open(outR1, 'w') as outfile1:
        with open(outR2, 'w') as outfile2:
            for read1, read2 in zip(file1, file2):
                Total += 1
                #look for valid barcode in forward read
                if len(forbarcodes) > 0:
                    BC, BCLabel = AlignBarcode(read1[1], forbarcodes, barcode_mismatch)
                    if BC == '':
                        NoBarcode += 1
                        continue
                if len(samples) > 0: #sample dictionary so enforce primers and barcodes from here
                    FwdPrimer = samples[BCLabel]['ForPrimer']
                    RevPrimer = samples[BCLabel]['RevPrimer']
                    foralign = edlib.align(FwdPrimer, read1[1], mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                    if foralign['editDistance'] < 0: #not found
                        NoPrimer += 1
                        continue                    
                    stringent = {}
                    stringent[BCLabel] = samples[BCLabel]['RevBarcode']
                    revBC, revBCLabel = AlignBarcode(read2[1], stringent, barcode_mismatch)
                    if revBC == '':
                        NoRevBarcode += 1
                        continue
                    #look for reverse primer in reverse read
                    revalign = edlib.align(RevPrimer, read2[1], mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                    if revalign['editDistance'] < 0: #not found
                        NoRevPrimer += 1
                        continue
                else:
                    #look for forward primer
                    foralign = edlib.align(fwdprimer, read1[1], mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                    if foralign['editDistance'] < 0: #not found
                        NoPrimer += 1
                        continue
                    if len(revbarcodes) > 0:
                        #look for valid revbarcodes
                        if len(revbarcodes) > 0:
                            revBC, revBCLabel = AlignBarcode(read2[1], revbarcodes, barcode_mismatch)
                            if revBC == '':
                                NoRevBarcode += 1
                                continue
                    #look for reverse primer in reverse read
                    revalign = edlib.align(revprimer, read2[1], mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
                    if revalign['editDistance'] < 0: #not found
                        NoRevPrimer += 1
                        continue
                #if get here, then all is well, construct new header and trim reads
                if BCLabel == revBCLabel:
                    label = BCLabel
                else:
                    label = BCLabel+':-:'+revBCLabel
                ForTrim = foralign["locations"][0][1]+1
                RevTrim = revalign["locations"][0][1]+1
                header = 'R_{:};barcodelabel={:};'.format(counter,label)
                outfile1.write('@%s\n%s\n+\n%s\n' % (header, read1[1][ForTrim:], read1[2][ForTrim:]))
                outfile2.write('@%s\n%s\n+\n%s\n' % (header, read2[1][RevTrim:], read2[2][RevTrim:]))
                counter += 1
                ValidSeqs += 1
    with open(stats, 'w') as statsout:
        statsout.write('%i,%i,%i,%i,%i,%i\n' % (Total, NoBarcode, NoPrimer, NoRevBarcode, NoRevPrimer, ValidSeqs)) 

def demuxIlluminaPE2(R1, R2, fwdprimer, revprimer, samples, forbarcodes, revbarcodes, barcode_mismatch, primer_mismatch, outR1, outR2, stats):
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest
    #function to loop through PE reads, renaming according to index
    file1 = FastqGeneralIterator(gzopen(R1))
    file2 = FastqGeneralIterator(gzopen(R2))
    counter = 1
    Total = 0
    NoBarcode = 0
    NoRevBarcode = 0
    NoPrimer = 0
    NoRevPrimer = 0
    ValidSeqs = 0
    with open(outR1, 'w') as outfile1:
        with open(outR2, 'w') as outfile2:
            for read1, read2 in zip(file1, file2):
                Total += 1
                #look for forward primer first, should all have primer and in correct orientation
                R1ForTrim = trimForPrimer(fwdprimer, read1[1], primer_mismatch)
                if R1ForTrim == 0: #no primer found
                    continue
                if len(forbarcodes) > 0: #search for barcode match in seq upstream of primer
                    R1BCtrim = R1ForTrim - len(fwdprimer)
                    BC, BCLabel = AlignBarcode2(read1[1][:R1BCtrim], forbarcodes, barcode_mismatch)
                    if BC == '':
                        NoBarcode += 1
                        continue                
                if len(samples) > 0: #sample dictionary so enforce primers and barcodes from here
                    FwdPrimer = samples[BCLabel]['ForPrimer']
                    RevPrimer = samples[BCLabel]['RevPrimer']             
                    stringent = {}
                    stringent[BCLabel] = samples[BCLabel]['RevBarcode']
                    #find rev primer in reverse read
                    R2ForTrim = trimForPrimer(RevPrimer, read2[1], primer_mismatch)
                    if R2ForTrim == 0:
                        continue
                    #look for reverse barcode
                    R2BCTrim = R2ForTrim - len(RevPrimer)
                    revBC, revBCLabel = AlignBarcode2(read2[1][:R2BCTrim], stringent, barcode_mismatch)
                    if revBC == '':
                        NoRevBarcode += 1
                        continue
                    #okay, found both primers and barcodes, now 1 more cleanup step trip revcomped primers
                    R1RevTrim = trimRevPrimer(RevPrimer, read1[1], primer_mismatch) 
                    R2RevTrim = trimRevPrimer(FwdPrimer, read2[1], primer_mismatch)
                else:
                    #no samples dictionary, so allow all combinations of matches
                    R2ForTrim = trimForPrimer(revprimer, read2[1], primer_mismatch)
                    if R2ForTrim == 0:
                        continue
                    if len(revbarcodes) > 0:
                        #look for reverse barcode
                        R2BCTrim = R2ForTrim - len(revprimer)
                        revBC, revBCLabel = AlignBarcode2(read2[1][:R2BCTrim], revbarcodes, barcode_mismatch)
                        if revBC == '':
                            NoRevBarcode += 1
                            continue
                    #okay, found both primers and barcodes, now 1 more cleanup step trip revcomped primers
                    R1RevTrim = trimRevPrimer(revprimer, read1[1], primer_mismatch) 
                    R2RevTrim = trimRevPrimer(fwdprimer, read2[1], primer_mismatch)                     

                #if get here, then all is well, construct new header and trim reads
                if BCLabel == revBCLabel:
                    label = BCLabel
                else:
                    label = BCLabel+':-:'+revBCLabel
                header = 'R_{:};barcodelabel={:};'.format(counter,label)
                outfile1.write('@%s\n%s\n+\n%s\n' % (header, read1[1][R1ForTrim:R1RevTrim], read1[2][R1ForTrim:R1RevTrim]))
                outfile2.write('@%s\n%s\n+\n%s\n' % (header, read2[1][R2ForTrim:R2RevTrim], read2[2][R2ForTrim:R2RevTrim]))
                counter += 1
                ValidSeqs += 1
    with open(stats, 'w') as statsout:
        statsout.write('%i,%i,%i,%i,%i,%i\n' % (Total, NoBarcode, NoPrimer, NoRevBarcode, NoRevPrimer, ValidSeqs)) 


def trimForPrimer(primer, seq, primer_mismatch):
    foralign = edlib.align(primer, seq, mode="HW", k=primer_mismatch, additionalEqualities=degenNuc)
    if foralign['editDistance'] < 0:
        return 0
    else:
        CutPos = foralign["locations"][0][1]+1
        return CutPos

def trimRevPrimer(primer, seq, primer_mismatch):
    revprimer = RevComp(primer)
    revalign = edlib.align(revprimer, seq, mode="HW", k=primer_mismatch, task="locations", additionalEqualities=degenNuc)
    if revalign['editDistance'] < 0:
        return len(seq)
    else:
        CutPos = revalign["locations"][0][0]
        return CutPos
        
def losslessTrim(input, fwdprimer, revprimer, mismatch, trimLen, padding, minlength, output):
    '''
    function to trim primers if found from SE reads
    and then trim/pad to a set length
    '''
    with open(output, 'w') as outfile:
        for title, seq, qual in FastqGeneralIterator(gzopen(input)):
            #sometimes primers sneek through the PE merging pipeline, check quickly again trim if found
            ForTrim = trimForPrimer(fwdprimer, seq, mismatch)
            RevTrim = trimRevPrimer(revprimer, seq, mismatch)
            Seq = seq[ForTrim:RevTrim]
            Qual = qual[ForTrim:RevTrim]
            if len(Seq) < int(minlength): #need this check here or primer dimers will get through
                continue
            if len(Seq) < int(trimLen) and padding == 'on':
                pad = int(trimLen) - len(Seq)
                SeqF = Seq + pad*'N'
                QualF = Qual + pad*'I'
            else:
                SeqF = Seq[:trimLen]
                QualF = Qual[:trimLen]
            outfile.write('@%s\n%s\n+\n%s\n' % (title, SeqF, QualF))


def checkBCinHeader(input):
    #read first header
    for title, seq, qual in FastqGeneralIterator(open(input)):
        header = title.split(' ')
        info = header[-1]
        if info.split(':')[-1].isdigit():
            return False
        else:
            return True
        break

def illuminaBCmismatch(R1, R2, index):
    remove = []
    for file in [R1,R2]:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            ID = title.split(' ')[0]
            BC = title.split(':')[-1]
            if BC != index:
                remove.append(ID)
    remove = set(remove)
    return remove

def fasta2barcodes(input, revcomp):
    BC = {}
    with open(input, 'r') as infile:
        for rec in SeqIO.parse(infile, 'fasta'):
            if not rec.id in BC:
                if revcomp:
                    Seq = RevComp(str(rec.seq))
                else:
                    Seq = str(rec.seq)
                BC[rec.id] = Seq
    return BC
            
def AlignBarcode(Seq, BarcodeDict, mismatch):
    besthit = []
    for BL in list(BarcodeDict.keys()):
        B = BarcodeDict[BL]
        #apparently people use N's in barcode sequences, this doesn't work well with
        #edlib, so if N then trim and then align, need to trim the seq as well
        if B.startswith('N'):
            origLen = len(B)
            B = B.lstrip('N')
            newLen = len(B)
            lenDiff = origLen - newLen
            newSeq = Seq[lenDiff:]
        else:
            newSeq = Seq
        align = edlib.align(B, newSeq, mode="SHW", k=mismatch, additionalEqualities=degenNuc)
        if align["editDistance"] < 0:
            continue
        elif align["editDistance"] == 0:
            return B, BL
        elif align["editDistance"] > 0 and mismatch == 0:
            continue
        else:
            if len(besthit) < 3:
                besthit = [B, BL, align["editDistance"]]
            else:
                if align["editDistance"] < int(besthit[2]):
                    besthit = [B, BL, align["editDistance"]]
    if len(besthit) == 3:
        return besthit[0], besthit[1]
    else:
        return "", ""
    
def AlignRevBarcode(Seq, BarcodeDict, mismatch):
    besthit = []
    for BL in list(BarcodeDict.keys()):
        B = BarcodeDict[BL]
        if B.endswith('N'):
            B = B.rstrip('N')
        align = edlib.align(B, Seq, mode="HW", k=mismatch, additionalEqualities=degenNuc)
        if align["editDistance"] < 0:
            continue
        elif align["editDistance"] == 0:
            return B, BL
        elif align["editDistance"] > 0 and mismatch == 0:
            continue
        else:
            if len(besthit) < 3:
                besthit = [B, BL, align["editDistance"]]
            else:
                if align["editDistance"] < int(besthit[2]):
                    besthit = [B, BL, align["editDistance"]]
    if len(besthit) == 3:
        return besthit[0], besthit[1]
    else:
        return "", ""

def AlignBarcode2(Seq, BarcodeDict, mismatch):
    besthit = []
    for BL in list(BarcodeDict.keys()):
        B = BarcodeDict[BL]
        #apparently people use N's in barcode sequences, this doesn't work well with
        #edlib, so if N then trim and then align, need to trim the seq as well
        if B.startswith('N'):
            origLen = len(B)
            B = B.lstrip('N')
            newLen = len(B)
            lenDiff = origLen - newLen
            newSeq = Seq[lenDiff:]
        else:
            newSeq = Seq
        align = edlib.align(B, newSeq, mode="HW", k=int(mismatch), additionalEqualities=degenNuc)
        if align["editDistance"] < 0:
            continue
        elif align["editDistance"] == 0:
            return B, BL
        elif align["editDistance"] > 0 and mismatch == 0:
            continue
        else:
            if len(besthit) < 3:
                besthit = [B, BL, align["editDistance"]]
            else:
                if align["editDistance"] < int(besthit[2]):
                    besthit = [B, BL, align["editDistance"]]
    if len(besthit) == 3:
        return besthit[0], besthit[1]
    else:
        return "", ""
    
def findFwdPrimer(primer, sequence, mismatch, equalities):
    #trim position
    TrimPos = None
    #search for match
    align = edlib.align(primer, sequence, mode="HW", k=int(mismatch), additionalEqualities=equalities)
    if align["editDistance"] >= 0: #we found a hit
        TrimPos = align["locations"][0][1]+1
    #return position will be None if not found
    return TrimPos

def findRevPrimer(primer, sequence, mismatch, equalities):
    #trim position
    TrimPos = None
    #search for match
    align = edlib.align(primer, sequence, mode="HW", task="locations", k=int(mismatch), additionalEqualities=equalities)
    if align["editDistance"] >= 0: #we found a hit
        TrimPos = align["locations"][0][0] 
    #return position will be None if not found
    return TrimPos
    
def MergeReadsSimple(R1, R2, tmpdir, outname, minlen, usearch, rescue, method):
    #check that num sequences is identical
    if not PEsanitycheck(R1, R2):
        log.error("%s and %s are not properly paired, exiting" % (R1, R2))
        sys.exit(1)
    #next run USEARCH/vsearch mergepe
    merge_out = os.path.join(tmpdir, outname + '.merged.fq')
    skip_for = os.path.join(tmpdir, outname + '.notmerged.R1.fq')
    report = os.path.join(tmpdir, outname +'.merge_report.txt')
    log.debug("Now merging PE reads")
    if method == 'usearch':
        cmd = [usearch, '-fastq_mergepairs', R1, '-reverse', R2, '-fastqout', merge_out, '-fastq_trunctail', '5', '-fastqout_notmerged_fwd', skip_for, '-minhsp', '12','-fastq_maxdiffs', '8', '-report', report, '-fastq_minmergelen', str(minlen), '-threads', '1']
    else:
        cmd = ['vsearch', '--fastq_mergepairs', R1, '--reverse', R2, '--fastqout', merge_out, '--fastqout_notmerged_fwd', skip_for, '--fastq_minmergelen', str(minlen), '--fastq_allowmergestagger', '--threads', '1']
    runSubprocess(cmd, log)
    #now concatenate files for downstream pre-process_illumina.py script
    final_out = os.path.join(tmpdir, outname)
    tmp_merge = os.path.join(tmpdir, outname+'.tmp')
    with open(tmp_merge, 'w') as cat_file:
        shutil.copyfileobj(open(merge_out,'r'), cat_file)
        if rescue == 'on':
            shutil.copyfileobj(open(skip_for,'r'), cat_file)
    #run phix removal
    #since most users have 32 bit usearch, check size of file, if > 3 GB, split into parts
    log.debug("Removing phix from %s" % outname)
    phixsize = getSize(tmp_merge)
    phixcount = countfastq(tmp_merge)
    log.debug('File Size: %i bytes' % phixsize)
    if phixsize > 3e9:
        log.debug('FASTQ > 3 GB, splitting FASTQ file into chunks to avoid potential memory problems with 32 bit usearch')
        phixdir = os.path.join(tmpdir, 'phix_'+str(os.getpid()))
        os.makedirs(phixdir)
        num = round(int((phixsize / 3e9))) + 1
        split_fastq(tmp_merge, phixcount, phixdir, int(num))
        for file in os.listdir(phixdir):
            if file.endswith(".fq"):
                output = os.path.join(phixdir, file+'.phix')
                file = os.path.join(phixdir, file)
                cmd = [usearch, '-filter_phix', file, '-output', output, '-threads', '1']
                runSubprocess(cmd, log)
        with open(final_out, 'w') as finalout:
            for file in os.listdir(phixdir):
                if file.endswith('.phix'):
                    with open(os.path.join(phixdir, file), 'r') as infile:
                        shutil.copyfileobj(infile, finalout)
        shutil.rmtree(phixdir)
    else:
        cmd = [usearch, '-filter_phix', tmp_merge, '-output', final_out, '-threads', '1']
        runSubprocess(cmd, log)
    #count output
    finalcount = countfastq(final_out)
    SafeRemove(merge_out)
    SafeRemove(skip_for)
    SafeRemove(tmp_merge)
    return phixcount, finalcount

    
def MergeReads(R1, R2, tmpdir, outname, read_length, minlen, usearch, rescue, method, index, mismatch):
    removelist = []
    if mismatch == 0 and index != '':
        if checkBCinHeader(R1):
            log.debug("Searching for index mismatches > 0: %s" % index) 
            removelist = illuminaBCmismatch(R1, R2, index)
            log.debug("Removing %i reads with index mismatch > 0" % len(removelist))  
    pretrim_R1 = os.path.join(tmpdir, outname + '.pretrim_R1.fq')
    pretrim_R2 = os.path.join(tmpdir, outname + '.pretrim_R2.fq')
    log.debug("Removing index 3prime bp 'A' from reads")
    trim3prime(R1, read_length, pretrim_R1, removelist)
    trim3prime(R2, read_length, pretrim_R2, removelist)  
    
    #check that num sequences is identical
    if not PEsanitycheck(pretrim_R1, pretrim_R2):
        log.error("%s and %s are not properly paired, exiting" % (R1, R2))
        sys.exit(1)
        
    #next run USEARCH/vsearch mergepe
    merge_out = os.path.join(tmpdir, outname + '.merged.fq')
    skip_for = os.path.join(tmpdir, outname + '.notmerged.R1.fq')
    report = os.path.join(tmpdir, outname +'.merge_report.txt')
    log.debug("Now merging PE reads")
    if method == 'usearch':
        cmd = [usearch, '-fastq_mergepairs', pretrim_R1, '-reverse', pretrim_R2, '-fastqout', merge_out, '-fastq_trunctail', '5', '-fastqout_notmerged_fwd', skip_for,'-minhsp', '12','-fastq_maxdiffs', '8', '-report', report, '-fastq_minmergelen', str(minlen)]
    else:
        cmd = ['vsearch', '--fastq_mergepairs', pretrim_R1, '--reverse', pretrim_R2, '--fastqout', merge_out, '--fastqout_notmerged_fwd', skip_for, '--fastq_minmergelen', str(minlen), '--fastq_allowmergestagger']
    runSubprocess(cmd, log)
    #now concatenate files for downstream pre-process_illumina.py script
    final_out = os.path.join(tmpdir, outname)
    tmp_merge = os.path.join(tmpdir, outname+'.tmp')
    with open(tmp_merge, 'w') as cat_file:
        shutil.copyfileobj(open(merge_out,'r'), cat_file)
        if rescue == 'on':
            shutil.copyfileobj(open(skip_for,'r'), cat_file)
    #run phix removal
    #since most users have 32 bit usearch, check size of file, if > 3 GB, split into parts
    log.debug("Removing phix from %s" % outname)
    phixsize = getSize(tmp_merge)
    phixcount = countfastq(tmp_merge)
    log.debug('File Size: %i bytes' % phixsize)
    if phixsize > 3e9:
        log.debug('FASTQ > 3 GB, splitting FASTQ file into chunks to avoid potential memory problems with 32 bit usearch')
        phixdir = os.path.join(tmpdir, 'phix_'+str(os.getpid()))
        os.makedirs(phixdir)
        num = round(int((phixsize / 3e9))) + 1
        split_fastq(tmp_merge, phixcount, phixdir, int(num))
        for file in os.listdir(phixdir):
            if file.endswith(".fq"):
                output = os.path.join(phixdir, file+'.phix')
                file = os.path.join(phixdir, file)
                cmd = [usearch, '-filter_phix', file, '-output', output]
                runSubprocess(cmd, log)
        with open(final_out, 'wb') as finalout:
            for file in os.listdir(phixdir):
                if file.endswith('.phix'):
                    with open(os.path.join(phixdir, file), 'r') as infile:
                        shutil.copyfileobj(infile, finalout)
        shutil.rmtree(phixdir)
    else:
        cmd = [usearch, '-filter_phix', tmp_merge, '-output', final_out]
        runSubprocess(cmd, log)
    #count output
    origcount = countfastq(R1)
    finalcount = countfastq(final_out)
    log.debug("Removed %i reads that were phiX" % (origcount - finalcount - len(removelist)))
    pct_out = finalcount / float(origcount) 
    #clean and close up intermediate files
    os.remove(merge_out)
    os.remove(pretrim_R1)
    os.remove(pretrim_R2)
    os.remove(skip_for)
    os.remove(tmp_merge)
    return log.info('{0:,}'.format(finalcount) + ' reads passed ('+'{0:.1%}'.format(pct_out)+')')

def validateorientation(tmp, reads, otus, output):
    orientcounts = os.path.join(tmp, 'orient.uc')
    cmd = ['vsearch', '--usearch_global', reads, '--db', otus, '--sizein', '--id', '0.97', '--strand', 'plus', '--uc', orientcounts]
    runSubprocess(cmd, log)
    OTUCounts = {}
    with open(orientcounts, 'r') as countdata:
        for line in countdata:
            line = line.rstrip()
            cols = line.split('\t')
            ID = cols[9]
            if ID == '*':
                continue
            size = cols[8].split('size=')[-1].replace(';', '')
            if not ID in OTUCounts:
                OTUCounts[ID] = int(size)
            else:
                OTUCounts[ID] += int(size)
    orientmap = os.path.join(tmp, 'orient-map.txt')
    cmd = ['vsearch', '--usearch_global', otus, '--db', otus, '--self', '--id', '0.95', '--strand', 'both', '--userout', orientmap, '--userfields', 'query+target+qstrand+id']
    runSubprocess(cmd, log)
    orient_remove = []
    keeper = []
    with open(orientmap, 'r') as selfmap:
        for line in selfmap:
            line = line.rstrip()
            cols = line.split('\t')
            if cols[2] == '-':
                qCount = OTUCounts.get(cols[0])
                tCount = OTUCounts.get(cols[1])
                if qCount > tCount:
                    if not cols[1] in orient_remove and not cols[1] in keeper:
                        orient_remove.append(cols[1])
                    if not cols[0] in keeper:
                        keeper.append(cols[0])
                else:
                    if not cols[0] in orient_remove and not cols[0]:
                        orient_remove.append(cols[0])
                    if not cols[1] in keeper:
                        keeper.append(cols[1])
    log.debug('Dropping {:,} OTUs: {:}'.format(len(orient_remove), ', '.join(orient_remove)))
    count = 0
    with open(output, 'w') as outfile:
        with open(otus, 'r') as infile:
            for rec in SeqIO.parse(infile, 'fasta'):
                if not rec.id in orient_remove:
                    count += 1
                    SeqIO.write(rec, outfile, 'fasta')
    return count, len(orient_remove)


def validateorientationDADA2(OTUCounts, otus, output):
    orientmap = 'orient-map.txt'
    cmd = ['vsearch', '--usearch_global', otus, '--db', otus, '--self', '--id', '0.95', '--strand', 'both', '--userout', orientmap, '--userfields', 'query+target+qstrand+id']
    runSubprocess(cmd, log)
    orient_remove = []
    keeper = []
    with open(orientmap, 'r') as selfmap:
        for line in selfmap:
            line = line.rstrip()
            cols = line.split('\t')
            if cols[2] == '-':
                qCount = OTUCounts.get(cols[0])
                tCount = OTUCounts.get(cols[1])
                if qCount > tCount:
                    if not cols[1] in orient_remove and not cols[1] in keeper:
                        orient_remove.append(cols[1])
                    if not cols[0] in keeper:
                        keeper.append(cols[0])
                else:
                    if not cols[0] in orient_remove and not cols[0]:
                        orient_remove.append(cols[0])
                    if not cols[1] in keeper:
                        keeper.append(cols[1])
    log.debug('Dropping {:,} OTUs: {:}'.format(len(orient_remove), ', '.join(natsorted(orient_remove))))
    count = 0
    with open(output, 'w') as outfile:
        with open(otus, 'r') as infile:
            for rec in SeqIO.parse(infile, 'fasta'):
                if not rec.id in orient_remove:
                    count += 1
                    SeqIO.write(rec, outfile, 'fasta')
    SafeRemove(orientmap)
    return count, len(orient_remove)


def dictFlip(input):
    #flip the list of dictionaries
    outDict = {}
    for k,v in input.items():
        for i in v:
            if not i in outDict:
                outDict[i] = k
            else:
                print("duplicate ID found: %s" % i)
    return outDict
    
def classifier2dict(input, pcutoff):
    ClassyDict = {}
    with open(input, 'r') as infile:
        for line in infile:
            cols = line.split('\t')
            ID = cols[0]
            tax = cols[1].split(',')
            passtax = []
            scores = []
            hit = False
            for i,level in enumerate(tax):
                if '(' in level:
                    score = level.split('(')[-1].replace(')', '')
                    if float(score) >= float(pcutoff):
                        hit = True
                        passtax.append(level.split('(')[0])
                        scores.append(score)
            if hit:
                if not ID in ClassyDict:
                    ClassyDict[ID] = (scores[-1], passtax)
    return ClassyDict

def usearchglobal2dict(input):
    GlobalDict = {}
    with open(input, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            cols = line.split('\t')
            ID = cols[0]
            if cols[1] != '*':
                tax = cols[1].split('tax=')[-1]
                tax = tax.split(',')
                pident = float(cols[-1]) / 100
                besthit = cols[1].split(';tax=')[0]
            else:
                tax = ['No hit']
                besthit = 'None'
                pident = 0.0000
            pident = "{0:.4f}".format(pident)
            #since we can have multiple hits with same pident, need to get best taxonomy
            if not ID in GlobalDict:
                GlobalDict[ID] = [(pident, besthit, '', tax,)]
            else:
                GlobalDict[ID].append((pident, besthit, '', tax))
    Results = {}
    for k,v in natsorted(GlobalDict.items()):
        mostTax = []
        lcaTax = []
        mt = 0
        for hit in v:
            if len(hit[-1]) > mt:
                mt = len(hit[-1])
        #loop through getting the ones with most tax info
        for x in v:
            if len(x[-1]) == mt:
                mostTax.append(x)
                lcaTax.append(x[-1])
        #now if mostTax has more than 1 hit, need to do LCA
        if len(mostTax) <= 1:
            Results[k] = mostTax[0]
        else:
            lcaResult = lcaTax[0]
            lcaFinal = []
            for x in lcaTax[1:]:
                s = set(x)
                lcaFinal = [z for z in lcaResult if z in s]
                lcaResult = lcaFinal
            if len(lcaResult) < mt:
                Results[k] = (mostTax[0][0], mostTax[0][1], 'LCA', lcaResult)
            else:
                Results[k] = (mostTax[0][0], mostTax[0][1], '',lcaResult)
    return Results

def bestclassifier(utax, sintax, otus):
    #both are dictionaries from classifier2dict, otus is a list of OTUs in order
    BestClassify = {}
    for otu in otus:
        sintaxhit, utaxhit = (None,)*2
        if otu in sintax:
            sintaxhit = sintax.get(otu) #this should be okay...
        if otu in utax:
            utaxhit = utax.get(otu) #returns tuple of (score, [taxlist]
        if sintaxhit and utaxhit: 
            if len(utaxhit[1]) > len(sintaxhit[1]):
                hit = (utaxhit[0], 'U', utaxhit[1])
            elif len(utaxhit[1]) == len(sintaxhit[1]):
                if float(utaxhit[0]) >= float(sintaxhit[0]):
                    hit = (utaxhit[0], 'U', utaxhit[1])
                else:
                    hit = (sintaxhit[0], 'S', sintaxhit[1])
            else:
                hit = (sintaxhit[0], 'S', sintaxhit[1])
        elif not sintaxhit and utaxhit:
            hit = (utaxhit[0], 'U', utaxhit[1])
        elif sintaxhit and not utaxhit:
            hit = (sintaxhit[0], 'S', sintaxhit[1])
        BestClassify[otu] = hit
    return BestClassify

def bestTaxonomy(usearch, classifier):
    Taxonomy = {}
    for k,v in natsorted(list(usearch.items())):
        globalTax = v[-1] #this should be a list
        classiTax = classifier.get(k)[-1] #also should be list
        pident = float(v[0])
        besthitID = v[1]
        LCA = v[2]
        discrep = 'S'
        if pident < 0.9700: #if global alignment hit is less than 97%, then default to classifier result
            method = classifier.get(k)[1].split()[0] #get first letter, either U or S
            score = float(classifier.get(k)[0])
            tax = ','.join(classiTax)
            fulltax = method+discrep+"|{0:.4f}".format(score)+'|'+besthitID+';'+tax     
        else: #should default to global alignment with option to update taxonomy from classifier if more information
            method = 'G'
            score = pident * 100
            if len(globalTax) < len(classiTax): #iterate through and make sure all levels match else stop where no longer matches
                tax = ','.join(classiTax)
                for lvl in globalTax:
                    if not lvl in classiTax: #now we have a problem
                        error = globalTax.index(lvl) #move one position backwards and keep only that level of taxonomy
                        discrep = 'D'
                        tax = ','.join(globalTax[:error])
                        break             
            else:
                tax = ','.join(globalTax)
            if discrep == 'S' and LCA == '':
                fulltax = method+discrep+"|{0:.1f}".format(score)+'|'+besthitID+';'+tax
            else:
                fulltax = method+discrep+"L|{0:.1f}".format(score)+'|'+besthitID+';'+tax
        Taxonomy[k] = fulltax
    return Taxonomy

def utax2qiime(input, output):
    domain = False
    with open(output, 'w') as outfile:
        outfile.write('#OTUID\ttaxonomy\n')
        with open(input, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                line = line.replace('\n', '')
                OTU = line.split('\t')[0]
                tax = line.split('\t')[1]
                try:
                    levels = tax.split(';')[1]
                except IndexError:
                    levels = '*'
                levels = levels.replace(',', ';')
                if levels.startswith('d:'):
                    domain = True
                    changes = ['d','p','c','o','f','g','s']
                else:
                    changes = ['k','p','c','o','f','g','s']
                for i in changes:
                    levels = levels.replace(i+':', i+'__')
                try:
                    levList = levels.split(';')
                except IndexError:
                    levList = [levels]
                #now loop through and add empty levels
                if not domain:
                    if not levList[0].startswith('k__'):
                        levList = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
                    if len(levList) < 2 and levList[0].startswith('k__'):
                        levList.extend(['p__', 'c__', 'o__', 'f__', 'g__', 's__'])
                    if len(levList) > 2 and not levList[2].startswith('c__'):
                        levList.insert(2,'c__')
                    if len(levList) > 3 and not levList[3].startswith('o__'):
                        levList.insert(3,'o__')
                    if len(levList) > 4 and not levList[4].startswith('f__'):
                        levList.insert(4,'f__')
                    if len(levList) > 5 and not levList[5].startswith('g__'):
                        levList.insert(5,'g__')
                    if len(levList) > 6 and not levList[6].startswith('s__'):
                        levList.insert(6,'s__')
                else:
                    if not levList[0].startswith('d__'):
                        levList = ['d__','p__', 'c__', 'o__', 'f__', 'g__', 's__']
                    if len(levList) < 2 and levList[0].startswith('d__'):
                        levList.extend(['p__', 'c__', 'o__', 'f__', 'g__', 's__'])                    
                    if len(levList) > 2 and not levList[2].startswith('c__'):
                        levList.insert(2,'c__')
                    if len(levList) > 3 and not levList[3].startswith('o__'):
                        levList.insert(3,'o__')
                    if len(levList) > 4 and not levList[4].startswith('f__'):
                        levList.insert(4,'f__')
                    if len(levList) > 5 and not levList[5].startswith('g__'):
                        levList.insert(5,'g__')
                    if len(levList) > 6 and not levList[6].startswith('s__'):
                        levList.insert(6,'s__')                             
                outfile.write('%s\t%s\n' % (OTU, ';'.join(levList)))

  
def barcodes2dict(input, mapDict, bcmismatch):
    #here expecting the index reads from illumina, create dictionary for naming?
    Results = {}
    NoMatch = []
    for title, seq, qual in FastqGeneralIterator(gzopen(input)):
        hit = [None, None, None]
        titlesplit = title.split(' ')
        readID = titlesplit[0]
        orient = titlesplit[1]
        if orient.startswith('2:'):
            seq = RevComp(seq)
        if seq in mapDict:
            BC = mapDict.get(seq)
            hit = [BC,seq, 0]
        else:
            for k,v in list(mapDict.items()):
                alignment = edlib.align(k, seq, mode="NW", k=bcmismatch)
                if alignment["editDistance"] < 0:
                    continue
                if hit[0]:
                    oldhit = hit[2]
                    if alignment["editDistance"] < oldhit:
                        hit = [v, k, alignment["editDistance"]]
                else:
                    hit = [v, k, alignment["editDistance"]]

            if not hit[0]: #no match, discard read
                NoMatch.append(readID)
                continue
        if not readID in Results:
            Results[readID] = (hit[0], hit[1], hit[2])
    return Results, NoMatch

def RevComp(s):
    rev_comp_lib = {'A':'T','C':'G','G':'C','T':'A','U':'A','M':'K','R':'Y','W':'W','S':'S','Y':'R','K':'M','V':'B','H':'D','D':'H','B':'V','X':'X','N':'N'}
    cseq = ''
    n = len(s)
    for i in range(0,n):
        c = s[n-i-1]
        cseq += rev_comp_lib[c.upper()]
    return cseq

def mapping2dict(input):
    #parse a qiime mapping file pull out seqs and ID into dictionary
    MapDict = {}
    IDs = []
    with open(input, 'r') as inputfile:
        for line in inputfile:
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            ID = cols[0]
            Seq = cols[1]
            if not Seq in MapDict:
                MapDict[Seq] = ID
            else:
                log.error("duplicate BC seq found %s: %s" % (Seq, ID))
            if not ID in IDs:
                IDs.append(ID)
            else:
                log.error("duplicate ID in mapping file: %s, exiting" (ID))
                sys.exit(1)
    return MapDict

def runMultiProgress(function, inputList, cpus, args=False):
    #setup pool
    p = multiprocessing.Pool(cpus)
    #setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, args=([i]), kwds={'args':args}))
    #refresh pbar every 5 seconds
    while True:
        incomplete_count = sum(1 for x in results if not x.ready())
        if incomplete_count == 0:
            break
        sys.stdout.write("     Progress: %.2f%% \r" % (float(tasks - incomplete_count) / tasks * 100))
        sys.stdout.flush()
        time.sleep(1)
    p.close()
    p.join()


def batch_iterator(iterator, batch_size):
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = next(iterator)
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
    if 'darwin' in sys.platform:
        stdoutformat = logging.Formatter(colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='[%b %d %I:%M %p]')      
    else:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='[%I:%M:%S %p]')
    fileformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='[%x %H:%M:%S]')
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
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    with open(output, 'w') as out:
        with open(input, 'r') as file:
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
    from Bio import SeqIO
    with open(input, 'r') as f:
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
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    seqs = {}
    with open(input, 'r') as file:
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
        with open(input, 'r') as fastq:
            for title, sequence, qual in FastqGeneralIterator(fastq):
                cols = title.split(';')
                header = 'R_'+str(count)+';'+cols[1]+';'
                count += 1
                out.write("@%s\n%s\n+\n%s\n" % (header, sequence, qual))

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

def fastarename(input, relabel, output):
    from Bio.SeqIO.FastaIO import FastaIterator
    with open(output, 'w') as outfile:
        counter = 1
        for record in FastaIterator(open(input)):
            newName = relabel+str(counter) 
            outfile.write(">%s\n%s\n" % (newName, record.seq))
            counter += 1

def fasta_strip_padding(file, output):
    from Bio.SeqIO.FastaIO import FastaIterator
    with open(output, 'w') as outputfile:
        for record in FastaIterator(gzopen(file)):
            Seq = record.seq.rstrip('N')   
            outputfile.write(">%s\n%s\n" % (record.id, Seq))

def fastq_strip_padding(file, output):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    with open(output, 'w') as outputfile:
        for title, seq, qual in FastqGeneralIterator(open(file)):
            Seq = seq.rstrip('N')
            Qual = qual[:len(Seq)]
            assert len(Seq) == len(Qual)    
            outputfile.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))
            
def ReverseComp(input, output):
    with open(output, 'w') as revcomp:
        with open(input, 'r') as fasta:
            for rec in SeqIO.parse(fasta, 'fasta'):
                revcomp.write(">%s\n%s\n" % (rec.id, rec.seq.reverse_complement()))

def guess_csv_dialect(header):
    """ completely arbitrary fn to detect the delimiter
    :type header: str
    :raise ValueError:
    :rtype: csv.Dialect
    """
    possible_delims = "\t,"
    lines = header.split("\n")
    if len(lines) < 2:
        raise ValueError("CSV header must contain at least 1 line")

    dialect = csv.Sniffer().sniff(header, delimiters=possible_delims)
    return dialect

    
def fasta2list(input):
    seqlist = []
    with open(input, 'r') as inseq:
        for rec in SeqIO.parse(inseq, 'fasta'):
            if not rec.description in seqlist:
                seqlist.append(rec.description)
    return seqlist

def updateMappingFile(mapfile, barcode_count, output):
    with open(output, 'w') as outfile:
        with open(mapfile, 'r') as infile:
            for line in infile:
                line = line.rstrip()
                cols = line.split('\t')
                if line.startswith('#Sample'): #header row, look for DemuxReads
                    if 'DemuxReads' in cols:
                        loc = cols.index('DemuxReads')
                    elif 'phinchID' in cols:
                        loc = cols.index('phinchID')+1
                        cols.insert(loc,'DemuxReads')
                    else:
                        cols.append('DemuxReads')
                        loc = cols.index('DemuxReads')
                    outfile.write('{:}\n'.format('\t'.join(cols)))
                else:
                    if cols[0] in barcode_count:
                        outfile.write('{:}\t{:}\t{:}\n'.format('\t'.join(cols[:loc]), str(barcode_count[cols[0]]), '\t'.join(cols[loc+1:])))
                        
def CreateGenericMappingFile(barcode_dict, revbarcode_dict, fwd_primer, rev_primer, output, barcodes_found):
    with open(output, 'w') as outfile:
        outfile.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tRevBarcodeSequence\tReversePrimer\tphinchID\tDemuxReads\tTreatment\n')
        for k,v in natsorted(barcodes_found.items()):
            sample = k
            Fsample,Rsample = (None,)*2
            if ':-:' in sample:
                Fsample, Rsample = sample.split(':-:')
            count = v
            forbarcode, revbarcode = ('no_data',)*2
            if Fsample:
                if Fsample in barcode_dict:
                    forbarcode = barcode_dict[Fsample]
            else:
                if sample in barcode_dict:
                    forbarcode = barcode_dict[sample]
            if Rsample:
                if Rsample in revbarcode_dict:
                    revbarcode = revbarcode_dict[Rsample]
            else:
                if sample in revbarcode_dict:
                    revbarcode = revbarcode_dict[sample]
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%i\t%s\n' % (sample, forbarcode, fwd_primer, revbarcode, rev_primer, sample, count, 'no_data'))
            

def CreateGenericMappingFileIllumina(samples, fwd_primer, rev_primer, output, barcodes):
    with open(output, 'w') as outfile:
        outfile.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tRevBarcodeSequence\tReversePrimer\tphinchID\tDemuxReads\tTreatment\n')
        for k,v in natsorted(list(samples.items())):
            count = barcodes.get(k, 0)
            if count > 0:
                if '-' in v:
                    forbarcode,revbarcode = v.split('-')
                else:
                    forbarcode = v
                    revbarcode = 'no_data'
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%i\t%s\n' % (k, forbarcode, fwd_primer, revbarcode, rev_primer, k, int(count), "no_data"))

def parseMappingFile(input, output):
    '''
    function to parse mapping file pull out primers and barcode sequences
    '''
    fwdprimer = ''
    revprimer = ''
    with open(output, 'w') as outfile:
        with open(input, 'r') as inputfile:
            for line in inputfile:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                outfile.write('>%s\n%s\n' % (cols[0], cols[1]))
                match = edlib.align(cols[1], cols[2], mode="HW", k=0)
                if match["editDistance"] == 0:
                    Trim = match["locations"][0][1]+1
                    if fwdprimer == '':
                        fwdprimer = cols[2][Trim:]
                        revprimer = cols[3]
    return (fwdprimer, revprimer)

def getMappingHeaderIndexes(input):
    IDx,FBCx,RBCx,FPx,RPx = (None,)*5
    with open(input, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('#'):
                cols = line.split('\t')
                if '#SampleID' in cols:
                    IDx = cols.index('#SampleID')
                if 'BarcodeSequence' in cols:
                    FBCx = cols.index('BarcodeSequence')
                if 'LinkerPrimerSequence' in cols:
                    FPx = cols.index('LinkerPrimerSequence')
                if 'RevBarcodeSequence' in cols:
                    RBCx = cols.index('RevBarcodeSequence')
                if 'ReversePrimer' in cols:
                    RPx = cols.index('ReversePrimer')
    return IDx, FBCx, FPx, RBCx, RPx

exampleMapFile='#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tRevBarcodeSequence\tReversePrimer\tphinchID\tTreatment\n'
                
def parseMappingFileNEW(input):
    '''
    function to parse mapping file pull out primers and barcode sequences
    '''
    results = {}
    ForBCDict = {}
    RevBCDict = {} 
    IDx, FBCx, FPx, RBCx, RPx = getMappingHeaderIndexes(input)
    if not any([IDx, FBCx, FPx, RPx]):
        log.error('Mapping file incorrectly formatted, headers should be (RevBarcodeSequence is optional):\n{:}'.format(exampleMapFile))
        sys.exit(1)
    if not RBCx:
        log.debug('Mapping file missing header: "RevBarcodeSequence", skipping reverse barcodes')
    with open(input, 'r') as inputfile:
        for line in inputfile:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            ID = cols[IDx]
            FBC = cols[FBCx]
            FP = cols[FPx]
            if FBC in FP: #barcode nested in primer_db
                loc = FP.index(FBC) + len(FBC)
                FP = FP[loc:]
            RP = cols[RPx]
            if RBCx:
                RBC = cols[RBCx]
                if RBC != '':
                    if RBC in RP:
                        loc = RP.index(RBC) + len(RBC)
                        RP = RP[loc:]
                else:
                    RBC = None
            else:
                RBC = None
            if not ID in results:
                results[ID] = {'ForBarcode': FBC, 'ForPrimer': FP, 'RevBarcode': RBC, 'RevPrimer': RP}
                ForBCDict[ID] = FBC
                if RBC:
                    RevBCDict[ID] = RBC
            else:
                log.error('Please fix duplicate SampleID detected in mapping file: {:}'.format(ID))
                sys.exit(1)
    return results, ForBCDict, RevBCDict, FP, RP

def parseMappingFileIllumina(input):
    fwdprimer = ''
    revprimer = ''
    samples = []
    with open(input, 'r') as inputfile:
        for line in inputfile:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if not cols[0] in samples:
                samples.append(cols[0])
            match = edlib.align(cols[1], cols[2], mode="HW", k=0)
            if match["editDistance"] == 0:
                Trim = match["locations"][0][1]+1
                if fwdprimer == '':
                    fwdprimer = cols[2][Trim:]
                    revprimer = cols[3]
            else:
                fwdprimer = cols[2]
                revprimer = cols[3]
        return (samples, fwdprimer, revprimer)
                          
def removefile(input):
    if os.path.isfile(input):
        os.remove(input)