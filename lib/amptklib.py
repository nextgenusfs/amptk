import sys, logging, csv, os, subprocess, multiprocessing, platform, time, shutil, inspect, gzip
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from natsort import natsorted
import revcomp_lib

ASCII = {'!':'0','"':'1','#':'2','$':'3','%':'4','&':'5',"'":'6','(':'7',')':'8','*':'9','+':'10',',':'11','-':'12','.':'13','/':'14','0':'15','1':'16','2':'17','3':'18','4':'19','5':'20','6':'21','7':'22','8':'23','9':'24',':':'25',';':'26','<':'27','=':'28','>':'29','?':'30','@':'31','A':'32','B':'33','C':'34','D':'35','E':'36','F':'37','G':'38','H':'39','I':'40','J':'41','K':'42','L':'43','M':'44','N':'45','O':'46','P':'47','Q':'48','R':'49','S':'50'}

primer_db = {'fITS7': 'GTGARTCATCGAATCTTTG', 'ITS4': 'TCCTCCGCTTATTGATATGC', 'ITS1-F': 'CTTGGTCATTTAGAGGAAGTAA', 'ITS2': 'GCTGCGTTCTTCATCGATGC', 'ITS3': 'GCATCGATGAAGAACGCAGC', 'ITS4-B': 'CAGGAGACTTGTACACGGTCCAG', 'ITS1': 'TCCGTAGGTGAACCTGCGG', 'LR0R': 'ACCCGCTGAACTTAAGC', 'LR2R': 'AAGAACTTTGAAAAGAG', 'JH-LS-369rc': 'CTTCCCTTTCAACAATTTCAC', '16S_V3': 'CCTACGGGNGGCWGCAG', '16S_V4': 'GACTACHVGGGTATCTAATCC', 'ITS3_KYO2': 'GATGAAGAACGYAGYRAA', 'COI-F': 'GGTCAACAAATCATAAAGATATTGG', 'COI-R': 'GGWACTAATCAATTTCCAAATCC', '515FB': 'GTGYCAGCMGCCGCGGTAA', '806RB': 'GGACTACNVGGGTWTCTAAT'}

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'
    
class gzopen(object):
   """Generic opener that decompresses gzipped files
   if needed. Encapsulates an open file or a GzipFile.
   Use the same way you would use 'open()'.
   """
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
   def next(self):
      return next(self.f)
      
def Funzip(input, output, cpus):
    '''
    function to unzip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '--decompress', '-c', '-p', str(cpus), input]
    elif which('bgzip'):
        cmd = ['bgzip', '--decompress', '-c', '-@', str(cpus), input]
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
    elif which('bgzip'):
        cmd = ['bgzip', '-c', '-@', str(cpus), input]
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
    elif which('bgzip'):
        cmd = ['bgzip', '-f', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '-f', input]
    try:
        runSubprocess(cmd, log)
    except NameError:
        subprocess.call(cmd)


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
    log.debug("Python Modules: %s" % ', '.join(results))

def fastalen2dict(input):
    Lengths = {}
    with open(input, 'rU') as infile:
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
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count
    
def countfastq(input):
    lines = sum(1 for line in gzopen(input))
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
    import pybam
    with open(output, 'w') as fastqout:
        with open(input, 'rb') as bamin:
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
    offsets = linepos[nstart:nstop]
    lines = []
    with gzopen(path) as inf:
        for offset in offsets:
            inf.seek(offset)
            lines.append(inf.readline())
    return lines     

def split_fastq(input, numseqs, outputdir, chunks):
    #get number of sequences and then number of sequences in each chunk
    numlines = numseqs*4
    n = numlines / chunks
    #make sure n is divisible by 4 (number of lines in fastq)
    if ( n % 4 ) != 0:
        n = ((n / 4) + 1) * 4
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
           
def MergeReads(R1, R2, tmpdir, outname, read_length, minlen, usearch, rescue, method, index, mismatch):
    removelist = []
    if mismatch == 0:
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
        cmd = ['vsearch', '--fastq_mergepairs', pretrim_R1, '--reverse', pretrim_R2, '--fastqout', merge_out, '--fastq_truncqual', '5', '--fastqout_notmerged_fwd', skip_for,'--fastq_maxdiffs', '8', '--fastq_minmergelen', str(minlen), '--fastq_allowmergestagger']
    runSubprocess(cmd, log)
    #now concatenate files for downstream pre-process_illumina.py script
    final_out = os.path.join(tmpdir, outname)
    tmp_merge = os.path.join(tmpdir, outname+'.tmp')
    with open(tmp_merge, 'w') as cat_file:
        shutil.copyfileobj(open(merge_out,'rU'), cat_file)
        if rescue == 'on':
            shutil.copyfileobj(open(skip_for,'rU'), cat_file)
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
                    with open(os.path.join(phixdir, file), 'rU') as infile:
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

def dictFlip(input):
    #flip the list of dictionaries
    outDict = {}
    for k,v in input.iteritems():
        for i in v:
            if not i in outDict:
                outDict[i] = k
            else:
                print "duplicate ID found: %s" % i
    return outDict
    
def fuzzymatch(seq1, seq2, num_errors):
    from Bio import pairwise2
    seq1_a, seq2_a, score, start, end = pairwise2.align.localms(seq1, seq2, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True, gap_char='-')[0]
    if start < 2: #align needs to start in first 2 bp
        if end < len(seq1)+2:
            match_region = seq2_a[start:end]
            seq_region = seq1_a[start:end]
            matches = sum((1 if s == match_region[i] else 0) for i, s in enumerate(seq_region))
            # too many errors -- no trimming
            if (len(seq1) - matches) <= int(num_errors):
                return (score, start, end)

def exactmatch(seq1, seq2):
    from Bio import pairwise2
    seq1_a, seq2_a, score, start, end = pairwise2.align.localms(seq1, seq2, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True, gap_char='-')[0]
    match_region = seq2_a[start:end]
    seq_region = seq1_a[start:end]
    matches = sum((1 if s == match_region[i] else 0) for i, s in enumerate(seq_region))
    # too many errors -- no trimming
    if (len(seq1) - matches) <= 0:
        return (score, start, end)
    else:
        return False

def classifier2dict(input, pcutoff):
    ClassyDict = {}
    with open(input, 'rU') as infile:
        for line in infile:
            cols = line.split('\t')
            ID = cols[0]
            tax = cols[1].split(',')
            passtax = []
            for level in tax:
                score = level.split('(')[-1].replace(')', '')
                if float(score) >= float(pcutoff):
                    passtax.append(level.split('(')[0])
                    oldscore = score
                else:
                    break
            ClassyDict[ID] = (oldscore, passtax)
    return ClassyDict

def usearchglobal2dict(input):
    GlobalDict = {}
    with open(input, 'rU') as infile:
        for line in infile:
            line = line.replace('\n', '')
            cols = line.split('\t')
            ID = cols[0]
            if cols[1] != '*':
                tax = cols[1].split('tax=')[-1]
                tax = tax.split(',')
                pident = float(cols[-1]) / 100
                besthit = cols[1].split(';tax=')[0]
            else:
                tax = 'No hit'
                besthit = 'None'
                pident = 0.0000
            pident = "{0:.4f}".format(pident)
            GlobalDict[ID] = (pident, besthit, tax)
    return GlobalDict

def bestclassifier(utax, sintax, otus):
    #both are dictionaries from classifier2dict, otus is a list of OTUs in order
    BestClassify = {}
    for otu in otus:
        utaxhit = utax.get(otu) #returns tuple of (score, [taxlist])
        sintaxhit = sintax.get(otu)
        if len(utaxhit[1]) > len(sintaxhit[1]):
            hit = (utaxhit[0], 'U', utaxhit[1])
        elif len(utaxhit[1]) == len(sintaxhit[1]):
            if float(utaxhit[0]) >= float(sintaxhit[0]):
                hit = (utaxhit[0], 'U', utaxhit[1])
            else:
                hit = (sintaxhit[0], 'S', sintaxhit[1])
        else:
            hit = (sintaxhit[0], 'S', sintaxhit[1])
        BestClassify[otu] = hit
    return BestClassify

def bestTaxonomy(usearch, classifier):
    Taxonomy = {}
    for k,v in natsorted(usearch.items()):
        globalTax = v[-1] #this should be a list
        classiTax = classifier.get(k)[-1] #also should be list
        pident = float(v[0])
        besthitID = v[1]
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
            fulltax = method+discrep+"|{0:.1f}".format(score)+'|'+besthitID+';'+tax
        Taxonomy[k] = fulltax
    return Taxonomy

def utax2qiime(input, output):
    with open(output, 'w') as outfile:
        outfile.write('#OTUID\ttaxonomy\n')
        with open(input, 'rU') as infile:
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
                changes = ['k','p','c','o','f','g','s']
                for i in changes:
                    levels = levels.replace(i+':', i+'__')
                try:
                    levList = levels.split(';')
                except IndexError:
                    levList = [levels]
                #now loop through and add empty levels
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
                outfile.write('%s\t%s\n' % (OTU, ';'.join(levList)))

  
def barcodes2dict(input, mapDict, bcmismatch):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    import edlib
    #here expecting the index reads from illumina, create dictionary for naming?
    Results = {}
    NoMatch = []
    for title, seq, qual in FastqGeneralIterator(gzopen(input)):
        hit = [None, None, None]
        titlesplit = title.split(' ')
        readID = titlesplit[0]
        orient = titlesplit[1]
        if orient.startswith('2:'):
            seq = revcomp_lib.RevComp(seq)
        if seq in mapDict:
            BC = mapDict.get(seq)
            hit = [BC,seq, 0]
        else:
            for k,v in mapDict.items():
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

def mapping2dict(input):
    #parse a qiime mapping file pull out seqs and ID into dictionary
    MapDict = {}
    IDs = []
    with open(input, 'rU') as inputfile:
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

def get_version():
    version = subprocess.Popen(['amptk', 'version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    return version

def get_usearch_version(usearch):
    try:
        version = subprocess.Popen([usearch, '-version'], stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    except OSError:
        log.error("%s not found in your PATH, exiting." % usearch)
        sys.exit(1)
    vers = version[0].replace('usearch v', '')
    vers2 = vers.split('_')[0]
    return vers2

def get_vsearch_version():
    version = subprocess.Popen(['vsearch', '--version'], stderr=subprocess.PIPE, stdout=subprocess.PIPE).communicate()
    for v in version:
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
    log.info("%s, USEARCH v%s, VSEARCH v%s" % (amptk_version, usearch_version, vsearch_version))

def checkRversion():
    #need to have R version > 3.2
    cmd = ['Rscript', '--vanilla', os.path.join(parentdir, 'util', 'check_version.R')]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    versions = stdout.replace(' \n', '')
    Rvers = versions.split(',')[0]
    dada2 = versions.split(',')[1]
    return (Rvers, dada2)

def runMultiProgress(function, inputList, cpus):
    #setup pool
    p = multiprocessing.Pool(cpus)
    #setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
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
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
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
    from Bio import SeqIO
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
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
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
        for title, seq, qual in FastqGeneralIterator(gzopen(file)):
            Seq = seq.rstrip('N')
            Qual = qual[:len(Seq)]
            assert len(Seq) == len(Qual)    
            outputfile.write("@%s\n%s\n+\n%s\n" % (title, Seq, Qual))
            
def ReverseComp(input, output):
    with open(output, 'w') as revcomp:
        with open(input, 'rU') as fasta:
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

def checkfastqsize(input):
    filesize = os.path.getsize(input)
    return filesize

def getCPUS():
    cores = multiprocessing.cpu_count()
    return cores
    
def fasta2list(input):
    seqlist = []
    with open(input, 'rU') as inseq:
        for rec in SeqIO.parse(inseq, 'fasta'):
            if not rec.description in seqlist:
                seqlist.append(rec.description)
    return seqlist
            

def CreateGenericMappingFile(barcode_fasta, fwd_primer, rev_primer, adapter, output, barcodes_found):
    mapDict = {}
    with open(barcode_fasta, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id in barcodes_found:
                if not rec.id in mapDict:
                    mapDict[rec.id] = (rec.seq, fwd_primer, rev_primer, rec.id, "no_data")
    with open(output, 'w') as outfile:
        outfile.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tphinchID\tTreatment\n')
        for k,v in natsorted(mapDict.items()):
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (k, v[0], adapter+v[0]+v[1], v[2], v[3], v[4]))

def CreateGenericMappingFileIllumina(samples, fwd_primer, rev_primer, output):
    with open(output, 'w') as outfile:
        outfile.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tphinchID\tTreatment\n')
        for k,v in natsorted(samples.items()):
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (k, v, fwd_primer, rev_primer, k, "no_data"))

def parseMappingFile(input, output):
    fwdprimer = ''
    revprimer = ''
    with open(output, 'w') as outfile:
        with open(input, 'rU') as inputfile:
            for line in inputfile:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                outfile.write('>%s\n%s\n' % (cols[0], cols[1]))
                match = exactmatch(cols[1], cols[2])
                if match:
                    if fwdprimer == '':
                        fwdprimer = cols[2][int(match[2]):]
                        revprimer = cols[3]
                else:
                    print "%s: barcode sequence not in LinkerPrimerSequence" % cols[0]
    return (fwdprimer, revprimer)

def parseMappingFileIllumina(input):
    fwdprimer = ''
    revprimer = ''
    samples = []
    with open(input, 'rU') as inputfile:
        for line in inputfile:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if not cols[0] in samples:
                samples.append(cols[0])
            match = exactmatch(cols[1], cols[2])
            if match:
                if fwdprimer == '':
                    fwdprimer = cols[2][int(match[2]):]
                    revprimer = cols[3]
            else:
                fwdprimer = cols[2]
                revprimer = cols[3]
        return (samples, fwdprimer, revprimer)
                          
def removefile(input):
    if os.path.isfile(input):
        os.remove(input)
        
        
'''        
def raup_crick(comm1, comm2, species_pool, reps=1000):
    s1, s2 = set(comm1), set(comm2)
    a1, a2 = len(comm1), len(comm2)
    ssobs = len(set.intersection(s1, s2))
    ssexp = [len(set.intersection(set(random.sample(species_pool, a1)),
                                  set(random.sample(species_pool, a2))
                                  )) 
             for _ in xrange(reps)]
    return ((len([s for s in ssexp if s > ssobs]) +
             len([s for s in ssexp if s == ssobs])/2.
             )/(len(ssexp)) - 0.5) * 2

 
def distance2mds(df, distance, type, output):
    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        from sklearn.metrics.pairwise import pairwise_distances
        from sklearn.manifold import MDS
        import seaborn as sns
        import matplotlib.pyplot as plt
    #run distance metric on matrix and then plot using NMDS
    num = len(df.index)
    data = np.array(df).astype(int)
    bc_dm = pairwise_distances(data, metric=distance)
    mds = MDS(n_components=2, metric=False, max_iter=999, dissimilarity='precomputed', n_init=10, verbose=0)
    result = mds.fit(bc_dm)
    coords = result.embedding_
    stress = 'stress=' + '{0:.4f}'.format(result.stress_)
    #get axis information and make square plus some padding
    xcoords = abs(maxabs(coords[:,0])) + 0.1
    ycoords = abs(maxabs(coords[:,1])) + 0.1
    #setup plot
    fig = plt.figure()
    #colors
    colorplot = sns.husl_palette(len(df), l=.5).as_hex()
    colorplot = [ str(x).upper() for x in colorplot ]
    for i in range(0,num):
        plt.plot(coords[i,0], coords[i,1], 'o', markersize=14, color=colorplot[i], label=df.index.values[i])
    plt.xlabel('NMDS axis 1')
    plt.ylabel('NMDS axis 2')
    plt.ylim(-ycoords,ycoords)
    plt.xlim(-xcoords,xcoords)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.title('NMDS analysis of '+type+' domains')
    plt.annotate(stress, xy=(1,0), xycoords='axes fraction', fontsize=12, ha='right', va='bottom')
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)
    
def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u

def _copy_array_if_base_present(a):
    """
    Copies the array if its base points to a parent array.
    """
    if a.base is not None:
        return a.copy()
    elif np.issubsctype(a, np.float32):
        return np.array(a, dtype=np.double)
    else:
        return a

def _convert_to_double(X):
    if X.dtype != np.double:
        X = X.astype(np.double)
    if not X.flags.contiguous:
        X = X.copy()
    return X


def raupcrick(u, v, n=1180, reps=1000):
    u = _validate_vector(u)
    v = _validate_vector(v)
    #length of vector is species pool?
    t = [1]*n
    #calculate alpha1, alpha2 (num observed in each site)
    a1, a2 = np.count_nonzero(u), np.count_nonzero(v)
    #calculate the number that are shared
    ssobs = np.sum(u == v)
    #now randomly sample the "species pool"
    idx = np.random.choice(np.arange(len(u)), size=reps)
    ssexp = np.sum(u[idx] == v[idx])
    
    return ((len([s for s in ssexp if s > ssobs]) +
             len([s for s in ssexp if s == ssobs])/2.
             )/(len(ssexp)) - 0.5) * 2 

def raupcrickdist(X, reps=1000):
    X = np.asarray(X, order='c')
    # The C code doesn't do striding.
    X = _copy_array_if_base_present(X)
    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')
    n, p = s
    D = np.zeros((n,n))
    #calc the species pool, proportion of sites each species occupies
    Z = []
    for i in X:
        Z.append(np.count_nonzero(i)/float(len(i)))
    Z = np.array(Z)
    #Loop through the ndarray
    for i in range(n):
        for j in range(n):
            s = 0
            for k in range(p):
                a1, a2 = np.count_nonzero(X[,p])
                s += (pts[i,k] - pts[j,k])**2
            m[i, j] = s**0.5
    return m
    for i in range(m): #looping through each row?
        for j in range(m):
            D[i] = 
    #calculate alpha1, alpha2 (num observed in each site)
    a1, a2 = np.count_nonzero(u), np.count_nonzero(v)
    #calculate the number that are shared
    ssobs = np.sum(u == v)
    #now randomly sample the "species pool"
    idx = np.random.choice(np.arange(len(u)), size=reps)
    ssexp = np.sum(u[idx] == v[idx])


def braycurt(u, v):
    u = _validate_vector(u)
    v = _validate_vector(v, dtype=np.float64)
    return abs(u - v).sum() / abs(u + v).sum()

def pdist(X, metric='euclidean', p=2, w=None, V=None, VI=None):
    # You can also call this as:
    #     Y = pdist(X, 'test_abc')
    # where 'abc' is the metric being tested.  This computes the distance
    # between all pairs of vectors in X using the distance metric 'abc' but with
    # a more succinct, verifiable, but less efficient implementation.

    X = np.asarray(X, order='c')

    # The C code doesn't do striding.
    X = _copy_array_if_base_present(X)

    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    dm = np.zeros((m * (m - 1)) // 2, dtype=np.double)

    wmink_names = ['wminkowski', 'wmi', 'wm', 'wpnorm']
    if w is None and (metric == wminkowski or metric in wmink_names):
        raise ValueError('weighted minkowski requires a weight '
                            'vector `w` to be given.')

    if callable(metric):
        if metric == minkowski:
            def dfun(u, v):
                return minkowski(u, v, p)
        elif metric == wminkowski:
            def dfun(u, v):
                return wminkowski(u, v, p, w)
        elif metric == seuclidean:
            def dfun(u, v):
                return seuclidean(u, v, V)
        elif metric == mahalanobis:
            def dfun(u, v):
                return mahalanobis(u, v, V)
        else:
            dfun = metric

        X = _convert_to_double(X)

        k = 0
        for i in xrange(0, m - 1):
            for j in xrange(i + 1, m):
                dm[k] = dfun(X[i], X[j])
                k = k + 1

    elif isinstance(metric, string_types):
        mstr = metric.lower()

        try:
            validate, pdist_fn = _SIMPLE_PDIST[mstr]
            X = validate(X)
            pdist_fn(X, dm)
            return dm
        except KeyError:
            pass

        if mstr in ['hamming', 'hamm', 'ha', 'h']:
            if X.dtype == bool:
                X = _convert_to_bool(X)
                _distance_wrap.pdist_hamming_bool_wrap(X, dm)
            else:
                X = _convert_to_double(X)
                _distance_wrap.pdist_hamming_wrap(X, dm)
        elif mstr in ['jaccard', 'jacc', 'ja', 'j']:
            if X.dtype == bool:
                X = _convert_to_bool(X)
                _distance_wrap.pdist_jaccard_bool_wrap(X, dm)
            else:
                X = _convert_to_double(X)
                _distance_wrap.pdist_jaccard_wrap(X, dm)
        elif mstr in ['minkowski', 'mi', 'm']:
            X = _convert_to_double(X)
            _distance_wrap.pdist_minkowski_wrap(X, dm, p)
        elif mstr in wmink_names:
            X = _convert_to_double(X)
            w = _convert_to_double(np.asarray(w))
            _distance_wrap.pdist_weighted_minkowski_wrap(X, dm, p, w)
        elif mstr in ['seuclidean', 'se', 's']:
            X = _convert_to_double(X)
            if V is not None:
                V = np.asarray(V, order='c')
                if V.dtype != np.double:
                    raise TypeError('Variance vector V must contain doubles.')
                if len(V.shape) != 1:
                    raise ValueError('Variance vector V must '
                                     'be one-dimensional.')
                if V.shape[0] != n:
                    raise ValueError('Variance vector V must be of the same '
                            'dimension as the vectors on which the distances '
                            'are computed.')
                # The C code doesn't do striding.
                VV = _copy_array_if_base_present(_convert_to_double(V))
            else:
                VV = np.var(X, axis=0, ddof=1)
            _distance_wrap.pdist_seuclidean_wrap(X, VV, dm)
        elif mstr in ['cosine', 'cos']:
            X = _convert_to_double(X)
            norms = _row_norms(X)
            _distance_wrap.pdist_cosine_wrap(X, dm, norms)
        elif mstr in ['old_cosine', 'old_cos']:
            X = _convert_to_double(X)
            norms = _row_norms(X)
            nV = norms.reshape(m, 1)
            # The numerator u * v
            nm = np.dot(X, X.T)
            # The denom. ||u||*||v||
            de = np.dot(nV, nV.T)
            dm = 1.0 - (nm / de)
            dm[xrange(0, m), xrange(0, m)] = 0.0
            dm = squareform(dm)
        elif mstr in ['correlation', 'co']:
            X = _convert_to_double(X)
            X2 = X - X.mean(1)[:, np.newaxis]
            norms = _row_norms(X2)
            _distance_wrap.pdist_cosine_wrap(X2, dm, norms)
        elif mstr in ['mahalanobis', 'mahal', 'mah']:
            X = _convert_to_double(X)
            if VI is not None:
                VI = _convert_to_double(np.asarray(VI, order='c'))
                VI = _copy_array_if_base_present(VI)
            else:
                if m <= n:
                    # There are fewer observations than the dimension of
                    # the observations.
                    raise ValueError("The number of observations (%d) is too "
                                     "small; the covariance matrix is "
                                     "singular. For observations with %d "
                                     "dimensions, at least %d observations "
                                     "are required." % (m, n, n + 1))
                V = np.atleast_2d(np.cov(X.T))
                VI = _convert_to_double(np.linalg.inv(V).T.copy())
            # (u-v)V^(-1)(u-v)^T
            _distance_wrap.pdist_mahalanobis_wrap(X, VI, dm)
        elif metric == 'test_euclidean':
            dm = pdist(X, euclidean)
        elif metric == 'test_sqeuclidean':
            if V is None:
                V = np.var(X, axis=0, ddof=1)
            else:
                V = np.asarray(V, order='c')
            dm = pdist(X, lambda u, v: seuclidean(u, v, V))
        elif metric == 'test_braycurtis':
            dm = pdist(X, braycurtis)
        elif metric == 'test_mahalanobis':
            if VI is None:
                V = np.cov(X.T)
                VI = np.linalg.inv(V)
            else:
                VI = np.asarray(VI, order='c')
            VI = _copy_array_if_base_present(VI)
            # (u-v)V^(-1)(u-v)^T
            dm = pdist(X, (lambda u, v: mahalanobis(u, v, VI)))
        elif metric == 'test_canberra':
            dm = pdist(X, canberra)
        elif metric == 'test_cityblock':
            dm = pdist(X, cityblock)
        elif metric == 'test_minkowski':
            dm = pdist(X, minkowski, p=p)
        elif metric == 'test_wminkowski':
            dm = pdist(X, wminkowski, p=p, w=w)
        elif metric == 'test_cosine':
            dm = pdist(X, cosine)
        elif metric == 'test_correlation':
            dm = pdist(X, correlation)
        elif metric == 'test_hamming':
            dm = pdist(X, hamming)
        elif metric == 'test_jaccard':
            dm = pdist(X, jaccard)
        elif metric == 'test_chebyshev' or metric == 'test_chebychev':
            dm = pdist(X, chebyshev)
        elif metric == 'test_yule':
            dm = pdist(X, yule)
        elif metric == 'test_matching':
            dm = pdist(X, matching)
        elif metric == 'test_dice':
            dm = pdist(X, dice)
        elif metric == 'test_kulsinski':
            dm = pdist(X, kulsinski)
        elif metric == 'test_rogerstanimoto':
            dm = pdist(X, rogerstanimoto)
        elif metric == 'test_russellrao':
            dm = pdist(X, russellrao)
        elif metric == 'test_sokalsneath':
            dm = pdist(X, sokalsneath)
        elif metric == 'test_sokalmichener':
            dm = pdist(X, sokalmichener)
        else:
            raise ValueError('Unknown Distance Metric: %s' % mstr)
    else:
        raise TypeError('2nd argument metric must be a string identifier '
                        'or a function.')
    return dm
'''