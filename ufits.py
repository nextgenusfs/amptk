#!/usr/bin/env python

#Wrapper script for UFITS package.

import sys, os, subprocess, inspect, tarfile, shutil, urllib2, urlparse
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList 
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

def fmtcols(mylist, cols):
    justify = []
    for i in range(0,cols):
        length = max(map(lambda x: len(x), mylist[i::cols]))
        length += 2
        ljust = map(lambda x: x.ljust(length), mylist[i::cols])
        justify.append(ljust)
    justify = flatten(justify)
    num_lines = len(mylist) / cols
    lines = (' '.join(justify[i::num_lines]) 
             for i in range(0,num_lines))
    return "\n".join(lines)

def download(url):
    file_name = "ufits_db.tar.gz"
    u = urllib2.urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
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


version = '0.4.1'

default_help = """
Usage:       ufits <command> <arguments>
version:     %s

Description: UFITS is a package of scripts to process fungal ITS amplicon data.  It uses the UPARSE algorithm for clustering
             and thus USEARCH8 is a dependency.
    
Process:     ion         pre-process Ion Torrent data (find barcodes, remove primers, trim/pad)
             illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, remove primers, trim/pad)
             illumina2   pre-process Illumina data from a single file (assumes Ion/454 read structure: <barcode><f_primer>READ)
             454         pre-process Roche 454 (pyrosequencing) data (find barcodes, remove primers, trim/pad)
             show        show number or reads per barcode from de-multiplexed data
             select      select reads from de-multiplexed data
             remove      remove reads from de-multiplexed data
             sample      sub-sample (rarify) de-multiplexed reads per sample
             
Clustering:  cluster     cluster OTUs (using UPARSE algorithm)
             cluster_ref closed/open reference based clustering
             filter      OTU table filtering
             taxonomy    Assign taxonomy to OTUs

Utilities:   summarize   Summarize Taxonomy (create OTU-like tables and/or stacked bar graphs for each level of taxonomy)
             funguild    Run FUNGuild (annotate OTUs with ecological information) 
             meta        pivot OTU table and append to meta data
             heatmap     Create heatmap from OTU table
             SRA         De-multiplex data and create meta data for NCBI SRA submission

Setup:       install     Download/install pre-formatted taxonomy DB (UNITE DB formatted for UFITS). Only need to run once.
             database    Format Reference Databases for Taxonomy
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
        """ % version

if len(sys.argv) > 1:
    if sys.argv[1] == 'ion':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script processes Ion Torrent PGM data for UFITS clustering.  The input to this script should be a 
             FASTQ file obtained from the Torrent Server analyzed with the `--disable-all-filters` flag to the 
             BaseCaller.  This script does the following: 1) finds Ion barcode sequences, 2) relabels headers with
             appropriate barcode name, 3) removes primer sequences, 4) trim/pad reads to a set length.
    
Arguments:   -i, --fastq,--bam   Input BAM or FASTQ file (Required)
             -o, --out           Output base name. Default: out
             -f, --fwd_primer    Forward primer sequence. Default: fITS7
             -r, --rev_primer    Reverse primer sequence Default: ITS4
             --barcode_fasta     FASTA file containing barcodes. Default: pgm_barcodes.fa
             -b, --barcodes      Barcodes used (list, e.g: 1,3,4,5,20). Default: all
             -n, --name_prefix   Prefix for re-naming reads. Default: R_
             -m, --min_len       Minimum length read to keep. Default: 50
             -l, --trim_len      Length to trim/pad reads. Default: 250
             --full_length       Keep only full length sequences.
             --primer_mismatch   Number of mismatches in primers to allow. Default: 2
             --cpus              Number of CPUs to use. Default: all
             --mult_samples      Combine multiple chip runs, name prefix for chip
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-process_ion.py')
            arguments.insert(0, cmd)
            arguments.append('--ion')
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'illumina2':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script takes Illumina MiSeq data that is not de-multiplexed and has read structure similar to Ion/454
             such that the reads are <barcode><fwd_primer>Read<rev_primer> for clustering using UFITS.  The default 
             behavior is to: 1) merge the PE reads using USEARCH, 2) find barcodes, 3)find and trim primers, 
             3) rename reads according to sample name, 4) trim/pad reads to a set length.
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             --reverse           Illumina PE reverse reads.
             -o, --out           Output base name. Default: out
             --barcode_fasta     FASTA file containing barcodes. Default: pgm_barcodes.fa
             -f, --fwd_primer    Forward primer sequence. Default: fITS7
             -r, --rev_primer    Reverse primer sequence Default: ITS4
             -n, --name_prefix   Prefix for re-naming reads. Default: R_
             -m, --min_len       Minimum length read to keep. Default: 50
             -l, --trim_len      Length to trim/pad reads. Default: 250
             --full_length       Keep only full length sequences.
             --cpus              Number of CPUs to use. Default: all
             -u, --usearch       USEARCH executable. Default: usearch8
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-process_ion.py')
            arguments.insert(0, cmd)
            arguments.append('--illumina')
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'illumina':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script takes a folder of Illumina MiSeq data that is already de-multiplexed and processes it for
             clustering using UFITS.  The default behavior is to: 1) merge the PE reads using USEARCH, 2) find and
             trim away primers, 3) rename reads according to sample name, 4) trim/pad reads to a set length.
    
Arguments:   -i, --fastq         Input folder of FASTQ files (Required)
             -o, --out           Output folder name. Default: ufits-data
             --reads             Paired-end or forward reads. Default: paired [paired, forward]
             --read_length       Illumina Read length (250 if 2 x 250 bp run). Default: 300 
             --rescue_forward    Rescue Forward Reads if PE do not merge, e.g. abnormally long amplicons
             -f, --fwd_primer    Forward primer sequence. Default: fITS7
             -r, --rev_primer    Reverse primer sequence Default: ITS4
             --require_primer    Require the Forward primer to be present. Default: on [on, off]
             -n, --name_prefix   Prefix for re-naming reads. Default: R_
             -m, --min_len       Minimum length read to keep. Default: 50
             -l, --trim_len      Length to trim/pad reads. Default: 250
             --full_length       Keep only full length sequences.
             --cpus              Number of CPUs to use. Default: all
             -u, --usearch       USEARCH executable. Default: usearch8
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-process_illumina_folder.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == '454':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script processes Roche 454 data for UFITS clustering.  The input to this script should be either a 
             SFF file, FASTA+QUAL files, or FASTQ file.  This script does the following: 1) finds barcode sequences, 
             2) relabels headers with appropriate barcode name, 3) removes primer sequences, 4) trim/pad reads to a set length.
    
Arguments:   -i, --sff, --fasta  Input file (SFF, FASTA, or FASTQ) (Required)
             -q, --qual          QUAL file (Required if -i is FASTA).
             -o, --out           Output base name. Default: out
             -f, --fwd_primer    Forward primer sequence. Default: fITS7
             -r, --rev_primer    Reverse primer sequence Default: ITS4
             --barcode_fasta     FASTA file containing barcodes. (Required)
             -n, --name_prefix   Prefix for re-naming reads. Default: R_
             -m, --min_len       Minimum length read to keep. Default: 50
             -l, --trim_len      Length to trim/pad reads. Default: 250
             --cpus              Number of CPUs to use. Default: all
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-process_ion.py')
            arguments.insert(0, cmd)
            arguments.append('--454')
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'cluster':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script is a "wrapper" for the UPARSE algorithm. FASTQ quality trimming via expected 
             errors and Dereplication are run in vsearch if installed otherwise defaults to Python 
             which allows for the use of datasets larger than 4GB.  
             Chimera filtering and UNOISE are also options.
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             -o, --out           Output base name. Default: out
             -e, --maxee         Expected error quality trimming. Default: 1.0
             -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
             -m, --minsize       Minimum size to keep (singleton filter). Default: 2
             --uchime_ref        Run Chimera filtering. Default: off [ITS1, ITS2, Full, 16S]
             --map_filtered      Map quality filtered reads back to OTUs. Default: off
             --unoise            Run De-noising pre-clustering (UNOISE). Default: off
             --size_annotations  Append size annotations to OTU names. Default: off
             -u, --usearch       USEARCH executable. Default: usearch8
             --cleanup           Remove intermediate files.
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-OTU_cluster.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)            
    elif sys.argv[1] == 'cluster_ref':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script first quality filters reads, dereplicates, and then runs chimera
             filtering.  OTUs are then picked via reference based clustering (closed)
             those that are > --id.  The rest of the data can then be clustered via
             de novo UPARSE and then reference clustered using UTAX.
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             -d, --db            Database Default: USEARCH
             --id                Percent ID for closed reference clustering. Default: 97
             --utax_db           UTAX formatted DB. Default: ITS2
             --utax_level        UTAX Taxonomy level to keep. Default: k [k,p,c,o,f,g,s]
             --utax_cutoff       UTAX confidence value threshold. Default: 0.8 [0 to 0.9]
             --mock              Mock community fasta file
             --closed_ref_only   Run only closed reference clustering.
             -o, --out           Output base name. Default: out
             -e, --maxee         Expected error quality trimming. Default: 1.0
             -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
             -m, --minsize       Minimum size to keep (singleton filter). Default: 2
             --uchime_ref        Run Chimera filtering. Default: off [ITS1, ITS2, Full, 16S]
             --map_filtered      Map quality filtered reads back to OTUs. Default: off
             -u, --usearch       USEARCH executable. Default: usearch8
             --cleanup           Remove intermediate files.
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-OTU_cluster_ref.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)            
    
    elif sys.argv[1] == 'filter':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script filters OTU table generated from the `ufits cluster` command and should be run on all datasets to combat
             barcode-switching or index-bleed (as high as 0.3 pct in MiSeq datasets, ~ 0.2 pct in Ion PGM datasets).  This script 
             works best when a spike-in control sequence is used, e.g. Synthetic Mock, although a mock is not required.
    
Required:    -i, --otu_table     OTU table
             -f, --fasta         OTU fasta
             
Optional:    -o, --out           Base name for output files. Default: use input basename
             -b, --mock_barcode  Name of barcode of mock community (Recommended)
             --mc                Mock community FASTA file. Default: ufits_synmock.fa 
             
Filtering    -n, --normalize     Normalize reads to number of reads per sample [y,n]. Default: y
             -p, --index_bleed   Filter index bleed between samples (percent). Default: 0.005
             -s, --subtract      Threshold to subtract from all OTUs (any number or auto). Default: 0
             -d, --delimiter     Delimiter of OTU tables. Default: csv  [csv, tsv]
             --min_reads_otu     Minimum number of reads for valid OTU from whole experiment. Default: 2
             --col_order         Column order (comma separated list). Default: sort naturally
             --keep_mock         Keep Spike-in mock community. Default: False
             -u, --usearch       USEARCH executable. Default: usearch8
             --show_stats        Show OTU stats on STDOUT   
             --cleanup           Remove intermediate files.
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-filter.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'select':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script filters de-multiplexed data (.demux.fq) to select only reads from samples provided
             in a text file, one name per line.
    
Required:    -i, --input     Input FASTQ file (.demux.fq)
             -l, --list      List of sample (barcode) names to keep, separate by space
             -f, --file      List of sample (barcode) names to keep in a file, one per line
             -o, --out       Output file name
             --format        File format for output file. Default: fastq [fastq, fasta]  
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'ufits-keep_samples.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)    
    elif sys.argv[1] == 'remove':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script filters de-multiplexed data (.demux.fq) to remove only reads from samples provided
             in a text file, one name per line.
    
Required:    -i, --input     Input FASTQ file (.demux.fq)
             -l, --list      List of sample (barcode) names to remove, separate by space
             -f, --file      List of sample (barcode) names to remove in a file, one per line
             -o, --out       Output file name
             --format        File format for output file. Default: fastq [fastq, fasta]
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'ufits-remove_samples.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)    
    elif sys.argv[1] == 'sample':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script sub-samples (rarifies) de-multiplexed data to equal number of reads per sample. For community
             analysis, this might not be appropriate as you are ignoring a portion of your data, however, there 
             might be some applications where it is useful.
    
Required:    -i, --input       Input FASTQ file
             -n, --num_reads   Number of reads to sub-sample to
             -o, --out         Output FASTQ file name      
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'ufits-barcode_rarify.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'meta':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script takes meta data file in CSV format (e.g. from excel) and an OTU table as input.  The first column
             of the meta data file must match the OTU table sample headers exactly.  It then pivots the OTU table and
             appends it to the meta data file.  
    
Required:    -i, --input       Input OTU table
             -m, --meta        Meta data table (csv format)
             -o, --out         Output (meta data + pivotted OTU table)
             --split_taxonomy  Make separate tables for groups of taxonomy [k,p,c,o,f,g]  
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'ufits-merge_metadata.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'show':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script takes de-multiplexed data (.demux.fq) as input and counts reads per barcode.
    
Required:    -i, --input     Input FASTQ file (.demux.fq)
             --quality_trim  Quality trim reads
             -e, --maxee     maxEE threshold for quality. Default: 1.0
             -l, --length    truncation length for trimming: Default: 250
             -o, --out       Output FASTQ file name (--quality_trim only)     
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'ufits-get_barcode_counts.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)    
    
    elif sys.argv[1] == 'funguild':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script takes OTU table as input and runs FUNGuild to assing functional annotation to an OTU
             based on the Guilds database.  Guilds script written by Zewei Song (2015).  
    
Options:     -i, --input        Input OTU table
             -d, --db           Database to use [fungi, nematode]. Default: fungi
             -o, --out          Output file basename.
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'Guilds_v1.0.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)    

    elif sys.argv[1] == 'heatmap':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script creates a heatmap from an OTU table.  Several settings are customizable.  Requires Matplotlib,
             numpy, and pandas.

Arguments:   -i, --otu_table     Input OTU table (Required)
             -o, --output        Output file (Required)
             --format            Image output format. Default: eps [eps, svg, png, pdf]
             --col_order         Column order (comma separated list). Default: sort naturally
             --font_size         Font Size for X/Y Axis labels. Default: 10
             --square            Maintain aspect ratio. Default: off
             --colors            Color Palette. Default: YlOrRd [many options, see matplotlib docs]
             --zero_color        Color for OTUs that are missing (zero). Default: white [white, lightgray, black, snow]
             --border_color      Color for border in-between cells. Default: black [black, white]
             --percent           Convert numbers to Percent of Sample. Default: off
             --min_num           Minimum reads per OTU to graph. Default: 1
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-heatmap.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)    
    elif sys.argv[1] == 'taxonomy':
        db_list = ['DB_name', 'DB_type', 'FASTA originated from', 'Fwd Primer', 'Rev Primer', 'Records']
        okay_list = []
        search_path = os.path.join(script_path, 'DB')
        for file in os.listdir(search_path):
            if file.endswith(".udb"):
                okay_list.append(file)
                info_file = file + '.txt'
                with open(os.path.join(search_path, info_file), 'rU') as info:
                    line = info.readlines()
                    line = [words for segments in line for words in segments.split()]
                    line.insert(0, file)
                    db_list.append(line)
        if len(db_list) < 7:
            db_print = "No DB configured, run 'ufits database' or 'ufits install' command."
        else:
            d = flatten(db_list)
            db_print = fmtcols(d, 6)
        
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script maps OTUs to taxonomy information and can append to an OTU table (optional).  By default the script
             uses a hybrid approach, e.g. gets taxonomy information from UTAX as well as global alignment hits from the larger
             UNITE-INSD database, and then parses both results to extract the most taxonomy information that it can at 
             'trustable' levels. UTAX results are used if BLAST-like search pct identity is less than 97 pct.  If pct identity
             is greater than 97 pct, the result with most taxonomy levels is retained.
    
Arguments:   -i, --otu_table     Input OTU table file (i.e. otu_table from ufits cluster) (Required)
             -f, --fasta         Input FASTA file (i.e. OTUs from ufits cluster) (Required)
             -o, --out           Base name for output file. Default: ufits-taxonomy.<method>.txt
             -m, --method        Taxonomy method. Default: hybrid [utax, usearch, hybrid, rdp, blast]
             --fasta_db          Alternative database of fasta sequenes to use for global alignment.
             --utax_db           UTAX formatted database. Default: ITS2.udb [See configured DB's below]
             --utax_cutoff       UTAX confidence value threshold. Default: 0.8 [0 to 0.9]
             --usearch_db        USEARCH formatted database. Default: USEARCH.udb
             --usearch_cutoff    USEARCH threshold percent identity. Default 0.7
             -r, --rdp           Path to RDP Classifier. Required if -m rdp
             --rdp_db            RDP Classifer DB set. [fungalits_unite, fungalits_warcup. fungallsu, 16srrna]  
             --rdp_cutoff        RDP Classifer confidence value threshold. Default: 0.8 [0 to 1.0]
             --local_blast       Local Blast database (full path) Default: NCBI remote nt database   
             --tax_filter        Remove OTUs from OTU table that do not match filter, i.e. Fungi to keep only fungi.
             -u, --usearch       USEARCH executable. Default: usearch8

Databases Configured: 
%s 
        """ % (sys.argv[1], version, db_print)

        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-assign_taxonomy.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)

        else:
            print help
            os._exit(1)                    
    elif sys.argv[1] == 'database':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Setup/Format reference database for ufits taxonomy command.
    
Arguments:   -i, --fasta         Input FASTA file (UNITE DB or UNITE+INSDC)
             -o, --out           Base Name for Output Files. Default: DB of ufits folder
             -f, --fwd_primer    Forward primer. Default: fITS7
             -r, --rev_primer    Reverse primer. Default: ITS4
             --format            Reformat FASTA headers to UTAX format. Default: unite2utax [unite2utax, rdp2utax, off]
             --drop_ns           Removal sequences that have > x N's. Default: 8
             --create_db         Create a DB. Default: usearch [utax, usearch]
             --skip_trimming     Keep full length sequences. Default: off
             --derep_fulllength  Remove identical sequences.
             --primer_mismatch   Max Primer Mismatch. Default: 4
             --keep_all          Keep Sequence if forward primer not found.
             --cpus              Number of CPUs to use. Default: all
             -u, --usearch       USEARCH executable. Default: usearch8       
        """ % (sys.argv[1], version)

        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-extract_region.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            try:
                outLocation = arguments.index('-o')
            except ValueError:
                outLocation = arguments.index('--out')
            outLocation = outLocation + 1
            arguments[outLocation] = os.path.join(script_path, 'DB', arguments[outLocation])
            subprocess.call(arguments)
        
        else:
            print help
            os._exit(1)   
    elif sys.argv[1] == 'summarize':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script traverses the taxonomy information and creates an OTU table for each
             level of taxonomy, i.e. Kingdom, Phylum, Class, etc.  Optionally, it will 
             create a Stacked Bar Graph for each taxonomy levels for each sample. Requires 
             Matplotlib, numpy, and pandas.
    
Arguments:   -i, --table     OTU Table containing Taxonomy information (Required)
             -o, --out       Base name for output files. Default: ufits-summary
             --graphs        Create stacked Bar Graphs.
             --format        Image output format. Default: eps [eps, svg, png, pdf]
             --percent       Convert numbers to Percent for Graphs. Default: off
             --font_size     Adjust font size for X-axis sample lables. Default: 8
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-summarize_taxonomy.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)        
        else:
            print help
            os._exit(1)  
    elif sys.argv[1] == 'install':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script downloads pre-formated UNITE and UNITE-INSD databases for use with the 
             `ufits taxonomy` command. By default UTAX is trained for ITS1, ITS2, and Full ITS. 
    
Arguments:   --install_unite     Install the UNITE Databases
             --force             Over-write existing databases
        """ % (sys.argv[1], version) 

        arguments = sys.argv[2:]
        if len(arguments) < 1:
            print help
            os._exit(1)
        else:
            if '--install_unite' in arguments:
                arguments.remove('--install_unite')
                #now check for UTAX and USEARCH in DB folder
                if os.path.isfile(os.path.join(script_path, 'DB', 'FULL.udb')):
                    if not '--force' in arguments:
                        print("A formated database was found, to overwrite use '--force'. You can add more custom databases by using the `ufits database` command.")
                        os._exit(1)                   
                #download database
                print "Downloading pre-formatted database"
                download('https://www.dropbox.com/s/mrrnqupzh43xyor/ufits_db.tar.gz?dl=1')
                tfile = tarfile.open("ufits_db.tar.gz", 'r:gz')
                tfile.extractall('ufits_db')
                for file in os.listdir('ufits_db'):
                    os.rename(os.path.join('ufits_db', file), os.path.join(script_path, 'DB', file))
                shutil.rmtree('ufits_db')
                os.remove('ufits_db.tar.gz')
                print "UFITS taxonomy database successfully installed"            

            else:
                print help
                os._exit(1)
    elif sys.argv[1] == 'SRA':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script aids in submitted your data to NCBI Sequence Read Archive (SRA) by splitting FASTQ file from Ion, 454, 
             or Illumina by barcode sequence into separate files for submission to SRA.  This ensures your data
             is minimally processed as only barcodes are removed.  Additionally, you can pass the --biosample argument
             with an NCBI biosample tab-delimited file and the script will auto-populate an SRA submission file.
    
Arguments:   -i, --input         Input FASTQ file or folder (Required)
             -o, --out           Output base name. Default: sra
             -b, --barcode_fasta Mulit-fasta file containing barcodes used.
             -s, --biosample     BioSample worksheet from NCBI (from confirmation email)
             -p, --platform      Sequencing platform. Defalt: ion (ion, illumina, 454)
             -f, --fwd_primer    Forward primer sequence. Default: fITS7
             -r, --rev_primer    Reverse primer sequence Default: ITS4
             -n, --names         CSV name mapping file, e.g. BC_1,NewName
             -d, --description   Paragraph description for SRA experimental design. Use quotes to wrap paragraph.
             --min_len           Minimum length read to keep after trimming barcodes. Default 50
             ---force            Overwrite directory with same name
        """ % (sys.argv[1], version)
   
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-fastq2sra.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
    elif sys.argv[1] == 'version':
        print "ufits v.%s" % version
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        os._exit(1)
    
    
else:
    print default_help
        