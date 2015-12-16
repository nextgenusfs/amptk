#!/usr/bin/env python

#Wrapper script for UFITS package.

import sys, os, subprocess, inspect
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


version = '0.2.10'

default_help = """
Usage:       ufits <command> <arguments>
version:     %s

Description: UFITS is a package of scripts to process fungal ITS amplicon data.  It uses the UPARSE algorithm for clustering
             and thus USEARCH8 is a dependency.
    
Command:     ion         pre-process Ion Torrent data (find barcodes, remove primers, trim/pad)
             illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, remove primers, trim/pad)
             illumina2   pre-process Illumina data from a single file (assumes Ion/454 read structure: <barcode><f_primer>READ)
             454         pre-process Roche 454 (pyrosequencing) data (find barcodes, remove primers, trim/pad)
             cluster     cluster OTUs (using UPARSE algorithm)
             filter      OTU table filtering
             taxonomy    Assign taxonomy to OTUs
             summarize   Summarize Taxonomy (create OTU-like tables and/or stacked bar graphs for each level of taxonomy)   
             heatmap     Create heatmap from OTU table

Setup:       install     Automated DB install (executes download and database commands for UNITE DBs). Only need to run once.
             download    Download Reference Databases
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
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             -o, --out           Output base name. Default: out
             -f, --fwd_primer    Forward primer sequence. Default: fITS7
             -r, --rev_primer    Reverse primer sequence Default: ITS4
             --barcode_fasta     FASTA file containing barcodes. Default: pgm_barcodes.fa
             -b, --barcodes      Barcodes used (list, e.g: 1,3,4,5,20). Default: all
             -n, --name_prefix   Prefix for re-naming reads. Default: R_
             -m, --min_len       Minimum length read to keep. Default: 50
             -l, --trim_len      Length to trim/pad reads. Default: 250
             --mult_samples      Combine multiple chip runs, name prefix for chip
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
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
             -u, --usearch       USEARCH executable. Default: usearch8
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
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
             --cpus              Number of CPUs to use. Default: all
             -u, --usearch       USEARCH executable. Default: usearch8
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
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
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
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

Description: Script is a "wrapper" for the UPARSE algorithm.  Modifications include support for a mock spike-in
             community.  FASTQ quality trimming via expected errors and Dereplication are run in Python which allows
             for the use of datasets larger than 4GB.  Chimera filtering and UNOISE are also options.
    
Arguments:   -i, --fastq         Input FASTQ file (Required)
             -o, --out           Output base name. Default: out
             -e, --maxee         Expected error quality trimming. Default: 1.0
             -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
             -m, --minsize       Minimum size to keep (singleton filter). Default: 2
             -l, --length        Length to trim reads. Default 250
             --mock              Name of spike-in mock community. Default: None
             --mc                Mock community FASTA file. Default: mock3
             --uchime_ref        Run Chimera filtering. Default: off [ITS1, ITS2, Full]
             --map_unfiltered    Map unfiltered reads back to OTUs. Default: off
             --unoise            Run De-noising pre-clustering (UNOISE). Default: off
             --size_annotations  Append size annotations to OTU names. Default: off
             -u, --usearch       USEARCH executable. Default: usearch8
             --cleanup           Remove intermediate files.
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
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
            
    elif sys.argv[1] == 'filter':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Script filters OTU table generated from the `ufits cluster` command.  There are two filters that
             are used: 1) index-bleed filter which combats barcode-switching or barcode mis-assignment and is
             controled by the --index_bleed argument, 2) the other filter is a subtraction filter meant to be
             used in conjuction with a spike-in mock community, it works asking the user to set a threshold based
             on the number of observed OTUs in the mock community versus expected number of OTUs.
    
Arguments:   -i, --otu_table     Input OTU table (Required)
             -b, --mock_barcode  Name of barcode of mock community (Required)
             -o, --out           Base name for output files. Default: use input basename
             -p, --index_bleed   Filter index bleed between samples (percent). Default: 0.1
             -d, --delimiter     Delimiter of OTU table. Default: tsv  [csv, tsv] 
             -n, --names         Change names of barcodes (CSV mapping file). Default: off
             --mc                Mock community FASTA file. Default: ufits_mock3.fa
             --col_order         Column order (comma separated list). Default: sort naturally
             --convert_binary    Convert OTU table to binary (1's and 0's). Default: off
             --trim_data         Filter the data. Default: on [on, off]
             --keep_mock         Keep Spike-in mock community. Default: False
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com      
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-mock_filter.py')
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
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
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
        #look in DB folder for databses
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
            db_print = "No DB configured, run 'ufits database' command for format database."
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
    
Arguments:   -i, --fasta         Input FASTA file (i.e. OTUs from ufits cluster) (Required)
             -o, --out           Base name for output file. Default: ufits-taxonomy.<method>.txt
             -m, --method        Taxonomy method. Default: hybrid [utax, usearch, hybrid, rdp, blast]
             --utax_db           UTAX formatted database. Default: ITS2.udb [See configured DB's below]
             --utax_cutoff       UTAX confidence value threshold. Default: 0.8 [0 to 0.9]
             --usearch_db        USEARCH formatted database. Default: USEARCH.udb
             --usearch_cutoff    USEARCH threshold percent identity. Default 0.7
             -r, --rdp           Path to RDP Classifier. Required if -m rdp
             --rdp_db            RDP Classifer DB set. [fungalits_unite, fungalits_warcup. fungallsu, 16srrna]  
             --rdp_cutoff        RDP Classifer confidence value threshold. Default: 0.8 [0 to 1.0]
             --local_blast       Local Blast database (full path) Default: NCBI remote nt database   
             --append_taxonomy   OTU table to append taxonomy. Default: none
             -u, --usearch       USEARCH executable. Default: usearch8

Databases Configured: 
%s
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
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
             --skip_trimming     Keep full length sequences. Default: off (not recommended)
             --primer_mismatch   Max Primer Mismatch. Default: 4
             --keep_all          Keep Sequence if forward primer not found.
             --cpus              Number of CPUs to use. Default: all
             -u, --usearch       USEARCH executable. Default: usearch8      
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
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
    
    elif sys.argv[1] == 'download':
        help = """
Usage:       ufits %s <arguments>
version:     %s

Description: Download reference databases for ufits taxonomy command.
    
Arguments:   -i, --input     Database to download. [unite, unite_insd, all]
     
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
        """ % (sys.argv[1], version)
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'ufits-download_db.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
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
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
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

Description: Script downloads the UNITE and UNITE-INSD databases and formats them for use with the 
             `ufits taxonomy` command. By default UTAX is trained for ITS1, ITS2, and Full ITS. 
    
Arguments:   --install_unite     Install the UNITE Databases 
             --primer_mismatch   Max Primer Mismatch. Default: 4
             -u, --usearch       USEARCH executable. Default: usearch8      
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
        """ % (sys.argv[1], version) 

        arguments = sys.argv[2:]
        if len(arguments) < 1:
            print help
            os._exit(1)
        else:
            if '--install_unite' in arguments:
                arguments.remove('--install_unite')
                primers = arguments
                #now check for UTAX and USEARCH in DB folder
                if os.path.isfile(os.path.join(script_path, 'DB', 'FULL.udb')):
                    print("A formated database was found, thus this command was already run.  You can add more custom databases by using the `ufits database` command.")
                    os._exit(1)
                #get UNITE general release
                arguments = ['-i', 'unite']
                cmd = os.path.join(script_path, 'bin', 'ufits-download_db.py')
                arguments.insert(0, cmd)
                exe = sys.executable
                arguments.insert(0, exe)
                subprocess.call(arguments)
                #now process DB for Full length and make UTAX db
                DB_location = os.path.join(script_path, 'DB', 'FULL')
                arguments1 = primers
                arguments1.append(['-i', 'sh_dynamic_01.08.2015.fasta', '-f', 'ITS1-F', '-r', 'ITS4', '--keep_all', '--derep_fulllength', '--create_db', 'utax', '-o', DB_location])
                arguments1 = flatten(arguments1)
                cmd1 = os.path.join(script_path, 'bin', 'ufits-extract_region.py')
                arguments1.insert(0, cmd1)
                arguments1.insert(0, exe)
                subprocess.call(arguments1)
                #now process DB for ITS1 and make UTAX db
                DB_location = os.path.join(script_path, 'DB', 'ITS1')
                arguments1 = primers
                arguments1.append(['-i', 'sh_dynamic_01.08.2015.fasta', '-f', 'ITS1-F', '-r', 'ITS2', '--keep_all', '--derep_fulllength', '--create_db', 'utax', '-o', DB_location])
                arguments1 = flatten(arguments1)
                cmd1 = os.path.join(script_path, 'bin', 'ufits-extract_region.py')
                arguments1.insert(0, cmd1)
                arguments1.insert(0, exe)
                subprocess.call(arguments1)
                #now process DB for ITS2 and make UTAX db
                DB_location = os.path.join(script_path, 'DB', 'ITS2')
                arguments1 = primers
                arguments1.append(['-i', 'sh_dynamic_01.08.2015.fasta', '-f', 'fITS7', '-r', 'ITS4','--derep_fulllength', '--create_db', 'utax', '-o', DB_location])
                arguments1 = flatten(arguments1)
                cmd1 = os.path.join(script_path, 'bin', 'ufits-extract_region.py')
                arguments1.insert(0, cmd1)
                arguments1.insert(0, exe)
                subprocess.call(arguments1)             
                #now get UNITE-INSD
                arguments = ['-i', 'unite_insd']
                cmd = os.path.join(script_path, 'bin', 'ufits-download_db.py')
                arguments.insert(0, cmd)
                exe = sys.executable
                arguments.insert(0, exe)
                subprocess.call(arguments)
                DB_location = os.path.join(script_path, 'DB', 'USEARCH')
                #now process DB
                arguments2 = primers
                arguments2.append(['-i', 'UNITE_public_01.08.2015.fasta', '--skip_trimming', '--derep_fulllength', '--create_db', 'usearch', '-o', DB_location])
                arguments2 = flatten(arguments2)
                cmd2 = os.path.join(script_path, 'bin', 'ufits-extract_region.py')
                arguments2.insert(0, cmd2)
                arguments2.insert(0, exe)
                subprocess.call(arguments2)
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
        