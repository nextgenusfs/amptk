#!/usr/bin/env python

#This script filters results from ufits-OTU_cluster.py
#written by Jon Palmer palmer.jona at gmail dot com

import os, argparse, inspect, subprocess, csv, logging, sys
from Bio import SeqIO
from natsort import natsorted

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(script_path)

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='ufits-mock_filter.py', usage="%(prog)s [options] -i otu_table.txt -b BC_27\n%(prog)s -h for help menu",
    description='''Script inspects output of ufits-OTU_cluster.py and 
    determines useful threshold for OTU output based on a spike-in 
    mock community.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--otu_table', required=True, help='Input OTU table')
parser.add_argument('-b','--mock_barcode', help='Barocde of Mock community')
parser.add_argument('-p','--index_bleed', dest="barcodebleed", default='0.1', help='Index Bleed filter')
parser.add_argument('-s','--subtract_threshold', default=0, type=int, help='Threshold to subtract')
parser.add_argument('--mc',default='mock3', help='Multi-FASTA mock community')
parser.add_argument('-d','--delimiter', default='tsv', choices=['csv','tsv'], help='Delimiter')
parser.add_argument('--col_order', dest="col_order", default="naturally", help='Provide comma separated list')
parser.add_argument('--convert_binary', dest="binary", action='store_true', help='Convert to binary')
parser.add_argument('--trim_data', default='on', choices=['on', 'off'], help='Threshold Trim Data')
parser.add_argument('-n', '--names', default='False', help='CSV mapping file BC,NewName')
parser.add_argument('--keep_mock', action='store_true', help='Keep mock sample in OTU table (Default: False)')
parser.add_argument('-o','--out', default='ufits-filter', help='Base output name')
args=parser.parse_args()

def try_int(x):
    try:
        return int(x)
    except ValueError:
        return x       
def try_subtract(x,y):
    try:
        z = x - y
        if z < 0:
            z = 0
        else:
            z = z
        return z
    except TypeError:
        return x       
def convert_binary(x):
    if isinstance(x, int):
        if x >= 1:
            return 1
        else:
            return 0
    else:
        return x
def greater_than(x):
    try:    
        if x >= blood_index:
            return x
        else:
            return 0
    except TypeError:
        return x

def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
    else:
        stdoutformat = logging.Formatter(colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='%b-%d-%Y %I:%M:%S %p')
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

#initialize script, log system info and usearch version
log.info("Operating system: %s" % sys.platform)
log.info("Loading OTU table: %s" % args.otu_table)
        
#check if otu_table is empty
check = os.stat(args.otu_table).st_size
if check == 0:
    log.error("Input OTU table is empty")
    os._exit(1)

#get base name of files
base = args.otu_table.split(".otu_table")
base = base[0]

#get default mock community value
if args.mc == "mock3":
    mock = os.path.join(parentdir, 'DB', 'ufits_mock3.fa')
elif args.mc == "mock2":
    mock = os.path.join(parentdir, 'DB', 'ufits_mock2.fa')
elif args.mc == "mock1":
    mock = os.path.join(parentdir, 'DB', 'ufits_mock1.fa')
else:
    mock = args.mc

#open mock community fasta and count records
mock_file = open(mock, "rU")
mock_ref_count = 0
for line in mock_file:
    if line.startswith (">"):
        mock_ref_count = mock_ref_count + 1
mock_file.close()

if not args.mock_barcode:
    log.info("No Mock Spike-in Barcode specified, will only apply index-bleed filter of %s percent" % args.barcodebleed)
    threshold = '0'
    sub_table = []
    keys = []
    trim_table = []
    binary_table = []
    re_orderTable = []
    finalTable = []
    line_count = 0
    if args.delimiter == 'csv':
        end = '.csv'
    if args.delimiter == 'tsv':
        end = '.txt'
    if args.out != 'ufits-filter':
        out_name = args.out + '.otu_table' + end
    else:
        if args.binary:
            out_name = base + '.filtered.otu_table.binary' + end
        else:
            out_name = base + '.filtered.otu_table' + end    
    file_out = open(out_name, "w")
    f2 = csv.reader(open(args.otu_table), delimiter='\t')
    for line in f2:
        line_count += 1
        sub_table.append([try_int(x) for x in line]) #convert to integers
    

    #shouldn't use subtract filter without mock community, but adding it here for another purpose
    tempTable = []
    if args.subtract_threshold != 0:
        log.info("Threshold subtracting table by %i" % args.subtract_threshold)
        for line in sub_table:
            tempTable.append([try_subtract(x,args.subtract_threshold) for x in line]) #subtract threshold 
    else:
        tempTable = sub_table
    #Now sort the table by row name naturally
    header = tempTable[:1] #pull out header row
    header = header[0] #just get first list, double check
        
     #convert header names (optional)
    if args.names != 'False':
        with open(args.names, 'rU') as infile:
            headRead = csv.reader(infile)
            headDict = {rows[0]:rows[1] for rows in headRead}
        header = [headDict[x] for x in header[1:]]
        header.insert(0,'OTUId')   
        mockBC = headDict.get(args.mock_barcode)
    else:
        mockBC = args.mock_barcode
    if args.col_order == "naturally":
        sortHead= natsorted(header)  #now sort the header row
    else:
        sortHead = args.col_order.split(",")
        sortHead.append('OTUId')

    otu = sortHead.index('OTUId') #pull out the OTU header
    sortHead.insert(0, sortHead.pop(otu)) #re-insert the OTU header in first column

    #sort the table without header row and then add back
    sortTable = natsorted(tempTable[1:]) #sort the table without header row
    sortTable.insert(0,header) #now add back header row

    #get index of the sorted header list 
    listIndex = []
    for i in sortHead:
        listIndex.append(header.index(i))

    #finally loop through the sorted table, and re-order the data
    count = 0
    for line in sortTable:
        count += 1
        lineList= []
        for item in listIndex:
            lineList.append(line[item])
        re_orderTable.append(lineList)        

    for line in re_orderTable:
        max_left = max(line[1:])
        if max_left >= 1:
            trim_table.append(line) #get rid of OTUs with only zeros
            keys.append(line[0])

    if args.barcodebleed != 'None':
        pct_bleed = float(args.barcodebleed) / 100
        temp_table = []
        for line in trim_table:
            if line[0] == 'OTUId':
                temp_table.append(line)
            else:
                subline = line[1:]
                OTU_sum = sum(subline)
                blood_index = OTU_sum * pct_bleed #set a threshold
                if OTU_sum > 0:
                    temp_table.append([greater_than(x) for x in line]) #filter blood_index
    else:
        temp_table = trim_table
    
    if args.binary:
        for line in temp_table:
            if line[0] == 'OTUId':
                binary_table.append(line)
            else:
                binary_table.append([convert_binary(x) for x in line]) #convert to binary
    else:
        binary_table = temp_table

    finalTable = binary_table

    #print out the final result
    for line in finalTable:
        if args.delimiter == 'tsv':
            join_delimiter = "\t"
        if args.delimiter == 'csv':
            join_delimiter = ","
        join_line = (join_delimiter.join(str(x) for x in line))
        file_out.write("%s\n" % join_line)
    file_out.close()

    #get new count of lines
    num_lines = sum(1 for line in open(out_name)) - 1
    line_count = line_count - 1
    #now lets write an updated OTU fasta file
    fasta_in = base + '.final.otus.fa'
    if args.out != 'ufits-filter':
        fasta_out = args.out + '.otus.fa'
    else:
        fasta_out = base + '.filtered.otus.fa'
    seqs_seen = []
    try:
        for record in SeqIO.parse(open(fasta_in, "rU"), "fasta"):
            if record.id in keys:
                seqs_seen.append(record)
        fasta_update = open(fasta_out, "w")
        SeqIO.write(seqs_seen, fasta_update, "fasta")
        fasta_update.close()
    except IOError:
        log.error( "Fasta file %s was not found, skipping writing fasta" % fasta_in)
    log.info("OTU table Filtering stats:\nOriginal OTUs: %i\nFiltered OTUs: %i" % (line_count, num_lines))

else:
    #load in OTU table, get only OTU column and mock, but also apply index bleed filter before subtraction filter, so print here
    KEEP_COLUMNS = ('OTUId', args.mock_barcode)
    f = csv.reader(open(args.otu_table), delimiter='\t')
    if args.barcodebleed != 'None':
        pct_bleed = float(args.barcodebleed) / 100
        int_table = []
        for line in f:
            int_table.append([try_int(x) for x in line]) #convert to integers
        temp_table = []
        for line in int_table:
            if line[0] == 'OTUId':
                temp_table.append(line)
            else:
                subline = line[1:]
                OTU_sum = sum(subline)
                blood_index = OTU_sum * pct_bleed #set a threshold
                if OTU_sum > 0:
                    temp_table.append([greater_than(x) for x in line]) #filter blood_index
    else:
        temp_table = f
    
    headers = None
    results = []
    for row in temp_table:
        if not headers:
            headers = []
            for i, col in enumerate(row):
                if col in KEEP_COLUMNS:
                    headers.append(i)
        else:
            results.append(tuple([row[i] for i in headers]))
    num_otus = 0
    mock_found = 0
    bad_otu = []
    good_otu = []
    for row in results:
        if int(row[1]) > 0:
            num_otus += 1
            if not "OTU" in row[0]:
                mock_found = mock_found + 1
                good_otu.append(int(row[1]))
            if "OTU" in row[0]:
                bad_otu.append(int(row[1]))
    spurious = num_otus - mock_found
    log.info("Summarizing data for %s" % (base))
    print "-------------------------------------------------------"
    print "Total OTUs in Mock:  %i" % (mock_ref_count)
    print "Total OTUs in %s:  %i" % (args.mock_barcode, num_otus)
    print "\nReal Mock OTUs found:  %i" % (mock_found)
    good_otu = sorted(good_otu, key=int)
    try:
        print "Range of counts from Real OTUs:  %i - %i" % (good_otu[-1], good_otu[0])
    except IndexError:
        print "\nThere does not appear to be Mock OTUs mapped in this table, run `ufits-OTU_cluster.py with the --mock parameter to generate a compatible OTU table.\n"
        os._exit(1)
    print "Lowest counts from Real OTUs:  %i, %i, %i" % (good_otu[0], good_otu[1], good_otu[2])
    print "\nSpurious OTUs found:  %i" % (spurious)
    if spurious != 0:
        bad_otu = sorted(bad_otu, key=int, reverse=True)
        print "Range of counts from Spurious OTUs:  %i - %i" % (bad_otu[0], bad_otu[-1])
        if spurious >= 3:
            print "Highest counts from Spurious OTUs:  %i, %i, %i" % (bad_otu[0], bad_otu[1], bad_otu[2])
        if spurious == 2:
            print "Highest counts from Spurious OTUs:  %i, %i" % (bad_otu[0], bad_otu[1])
        if spurious == 1:
            print "Highest count from Spurious OTUs:  %i" % (bad_otu[0])

    if args.trim_data == 'on':
            threshold = raw_input("\nEnter threshold value to trim data:  ")
            print "-------------------------------------------------------"
            num = int(threshold)
            new_table = []
            sub_table = []
            keys = []
            trim_table = []
            binary_table = []
            re_orderTable = []
            finalTable = []
            line_count = 0
            if args.delimiter == 'csv':
                end = '.csv'
            if args.delimiter == 'tsv':
                end = '.txt'
            if args.out != 'ufits-filter':
                out_name = args.out + '.otu_table' + end
            else:
                if args.binary:
                    out_name = base + '.filtered_' + threshold + '.otu_table.binary' + end
                else:
                    out_name = base + '.filtered_' + threshold + '.otu_table' + end
            
            file_out = open(out_name, "w")
            f2 = csv.reader(open(args.otu_table), delimiter='\t')
            for line in f2:
                line_count += 1
                new_table.append([try_int(x) for x in line]) #convert to integers
            
            if args.barcodebleed != 'None':
                pct_bleed = float(args.barcodebleed) / 100
                temp_table = []
                for line in new_table:
                    if line[0] == 'OTUId':
                        temp_table.append(line)
                    else:
                        subline = line[1:]
                        OTU_sum = sum(subline)
                        blood_index = OTU_sum * pct_bleed #set a threshold
                        if OTU_sum > 0:
                            temp_table.append([greater_than(x) for x in line]) #filter blood_index
            else:
                temp_table = new_table  
            
            for line in temp_table:
                sub_table.append([try_subtract(x,num) for x in line]) #subtract threshold          
        
            #Now sort the table by row name naturally
            header = sub_table[:1] #pull out header row
            header = header[0] #just get first list, double check
                
             #convert header names (optional)
            if args.names != 'False':
                with open(args.names, 'rU') as infile:
                    headRead = csv.reader(infile)
                    headDict = {rows[0]:rows[1] for rows in headRead}
                header = [headDict[x] for x in header[1:]]
                header.insert(0,'OTUId')   
                mockBC = headDict.get(args.mock_barcode)
            else:
                mockBC = args.mock_barcode
            if args.col_order == "naturally":
                sortHead= natsorted(header)  #now sort the header row
            else:
                sortHead = args.col_order.split(",")
                sortHead.append('OTUId')
        
            otu = sortHead.index('OTUId') #pull out the OTU header
            sortHead.insert(0, sortHead.pop(otu)) #re-insert the OTU header in first column
        
            #remove mock spike in control
            if not args.keep_mock:
                sortHead.remove(mockBC)
            
            #sort the table without header row and then add back
            sortTable = natsorted(sub_table[1:]) #sort the table without header row
            sortTable.insert(0,header) #now add back header row

            #get index of the sorted header list 
            listIndex = []
            for i in sortHead:
                listIndex.append(header.index(i))

            #finally loop through the sorted table, and re-order the data
            count = 0
            for line in sortTable:
                count += 1
                lineList= []
                for item in listIndex:
                    lineList.append(line[item])
                re_orderTable.append(lineList)        
        
            for line in re_orderTable:
                max_left = max(line[1:])
                if max_left >= 1:
                    trim_table.append(line) #get rid of OTUs with only zeros
                    keys.append(line[0])

            
            if args.binary:
                for line in trim_table:
                    if line[0] == 'OTUId':
                        binary_table.append(line)
                    else:
                        binary_table.append([convert_binary(x) for x in line]) #convert to binary
            else:
                binary_table = trim_table
        
            finalTable = binary_table
        
            #print out the final result
            for line in finalTable:
                if args.delimiter == 'tsv':
                    join_delimiter = "\t"
                if args.delimiter == 'csv':
                    join_delimiter = ","
                join_line = (join_delimiter.join(str(x) for x in line))
                file_out.write("%s\n" % join_line)
            file_out.close()
        
            #get new count of lines
            num_lines = sum(1 for line in open(out_name)) - 1
            line_count = line_count - 1
            #now lets write an updated OTU fasta file
            fasta_in = base + '.final.otus.fa'
            if args.out != 'ufits-filter':
                fasta_out = args.out + '.otus.fa'
            else:
                fasta_out = base + '.filtered_' + threshold + '.otus.fa'
            seqs_seen = []
            try:
                for record in SeqIO.parse(open(fasta_in, "rU"), "fasta"):
                    if record.id in keys:
                        seqs_seen.append(record)
                fasta_update = open(fasta_out, "w")
                SeqIO.write(seqs_seen, fasta_update, "fasta")
                fasta_update.close()
            except IOError:
                log.error("Fasta file %s was not found, skipping writing fasta" % fasta_in)
            log.info("OTU table Filtering stats:\nOTU table has been filtered to %i:\nOriginal OTUs: %i\nFiltered OTUs: %i" % (num, line_count, num_lines))
print "-------------------------------------------------------"
log.info("New OTU table created: %s" % out_name)
log.info("New FASTA OTU file: %s" % fasta_out)
print "-------------------------------------------------------"

if 'win32' in sys.platform:
    print "\nExample of next cmd: ufits taxonomy -i %s --append_taxonomy\n" % (fasta_out)
else:
    print colr.WARN + "\nExample of next cmd:" + colr.END +  " ufits taxonomy -i %s --append_taxonomy\n" % (fasta_out)
