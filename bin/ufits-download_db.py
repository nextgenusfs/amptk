#!/usr/bin/env python

import sys, os, re, argparse, logging, gzip, urllib2, shutil, urlparse, zipfile

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)

class col:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

parser=argparse.ArgumentParser(prog='ufits-download_db.py', usage="%(prog)s [options] -i <DB name>",
    description='''Script downloads ITS reference databases.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', dest='input', required=True, choices=['unite', 'unite_insd', 'rtl_ITS', 'all'], help='Download Reference Database (Required)')
args=parser.parse_args()

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


def download(url, fileName=None):
    global filesaved
    def getFileName(url,openUrl):
        if 'Content-Disposition' in openUrl.info():
            # If the response has Content-Disposition, try to get filename from it
            cd = dict(map(
                lambda x: x.strip().split('=') if '=' in x else (x.strip(),''),
                openUrl.info()['Content-Disposition'].split(';')))
            if 'filename' in cd:
                filename = cd['filename'].strip("\"'")
                if filename: return filename
        # if no filename was found above, parse it out of the final URL.
        return os.path.basename(urlparse.urlsplit(openUrl.url)[2])

    r = urllib2.urlopen(urllib2.Request(url))
    try:
        fileName = fileName or getFileName(url,r)
        with open(fileName, 'wb') as f:
            shutil.copyfileobj(r,f)
    finally:
        filesaved = fileName
        r.close()

def getUNITE(url):
    log.info("Downloading UNITE DB")
    download(url)
    log.info("Download finished, now unzipping: %s" % filesaved)
    with zipfile.ZipFile(filesaved, "r") as z:
        z.extractall()
    log.info("Cleaning up, file is saved as: %s" % unite_short)
    shutil.move(os.path.join('developer', unite_keep), unite_short)
    shutil.rmtree('developer')
    os.remove(filesaved)
    os.remove(unite_delete)
    
def getINSD(url):
    log.info("Downloading UNITE/INSDC DB")
    download(url)
    log.info("Download finished, now unzipping: %s" % filesaved)
    with zipfile.ZipFile(filesaved, "r") as z:
        z.extract(insd_keep)
    log.info("Cleaning up, file saved as: %s" % insd_keep)
    os.remove(filesaved)

def getRTL(url):
    log.info("Downloading NCBI RefSeq RTL ITS database")
    download(url)
    log.info("Download finished, file saved as: %s" % rtl_keep)

log_name = ('ufits_db.log')
if os.path.isfile(log_name):
    os.remove(log_name)

setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
log.debug(cmd_args)
print "-------------------------------------------------------"

#locations
unite_url = 'https://unite.ut.ee/sh_files/sh_general_release_s_01.08.2015.zip'
unite_keep = 'sh_general_release_dynamic_s_01.08.2015_dev.fasta'
unite_short = 'sh_dynamic_01.08.2015.fasta'
unite_delete = 'sh_general_release_dynamic_s_01.08.2015.fasta'
insd_url = 'https://unite.ut.ee/sh_files/UNITE_public_01.08.2015.fasta.zip'
insd_keep = 'UNITE_public_01.08.2015.fasta'
rtl_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/TARGET/ITS_rRNA/Fungi/Fungi.fna'
rtl_keep = 'Fungi.fna'

if args.input == 'unite':
    getUNITE(unite_url)

elif args.input == 'unite_insd':
    getINSD(insd_url)

elif args.input == 'rtl_ITS':
    getRTL(rtl_url)

elif args.input == 'all':
    getUNITE(unite_url)
    getINSD(insd_url)
    getRTL(rtl_url)
    
