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
    maxwidth = max(map(lambda x: len(x), mylist))
    justifyList = map(lambda x: x.ljust(maxwidth), mylist)
    lines = (' '.join(justifyList[i:i+cols]) 
             for i in xrange(0,len(justifyList),cols))
    return "\n".join(lines)

version = '0.2.5'

default_help = """
Usage:      ufits <command> <arguments>
version:    %s
    
Command:    ion         pre-process Ion Torrent data (find barcodes, remove primers, trim/pad)
            illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, remove primers, trim/pad)
            cluster     cluster OTUs (using UPARSE algorithm)
            filter      OTU table filtering
            taxonomy    Assign taxonomy to OTUs
            summarize   Summarize Taxonomy (create stacked bar graph and data tables)   
            heatmap     Create heatmap from OTU table

Setup:      download    Download Reference Databases
            database    Format Reference Databases for Taxonomy
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
        """ % version

if len(sys.argv) > 1:
    if sys.argv[1] == 'ion':
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --fastq         Input FASTQ file (Required)
            -o, --out           Output base name. Default: out
            -f, --fwd_primer    Forward primer sequence. Default: AGTGARTCATCGAATCTTTG (fITS7)
            -r, --rev_primer    Reverse primer sequence Default: TCCTCCGCTTATTGATATGC (ITS4)
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
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit
    elif sys.argv[1] == 'illumina':
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --fastq         Input FASTQ file (Required)
            -o, --out           Output folder name. Default: ufits-data
            --reads             Paired-end or forward reads. Default: paired [paired, forward]
            -f, --fwd_primer    Forward primer sequence. Default: GTGARTCATCGAATCTTTG (fITS7)
            -r, --rev_primer    Reverse primer sequence Default: TCCTCCGCTTATTGATATGC (ITS4)
            --require_primer    Require the Forward primer to be present. Default: on [on, off]
            -n, --name_prefix   Prefix for re-naming reads. Default: R_
            -m, --min_len       Minimum length read to keep. Default: 50
            -l, --trim_len      Length to trim/pad reads. Default: 250
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
            sys.exit
    
    elif sys.argv[1] == 'cluster':
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --fastq         Input FASTQ file (Required)
            -o, --out           Output base name. Default: out
            -e, --maxee         Expected error quality trimming. Default: 1.0
            -p, --pct_otu       OTU Clustering Radius (percent). Default: 97
            -m, --minsize       Minimum size to keep (singleton filter). Default: 2
            -l, --length        Length to trim reads. Default 250
            --mock              Name of spike-in mock community. Default: None
            --mc                Mock community FASTA file. Default: ufits_mock3.fa
            --uchime_ref        Run Chimera filtering. Default: off [ITS1, ITS2, Full]
            --map_unfiltered    Map unfiltered reads back to OTUs. Default: off
            --unoise            Run De-noising pre-clustering (UNOISE). Default: off
            --size_annotations  Append size annotations to OTU names. Default: off
            -u, --usearch       USEARCH executable. Default: usearch8
            
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
            sys.exit
            
    elif sys.argv[1] == 'filter':
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --otu_table     Input OTU table (Required)
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
            sys.exit
    
    elif sys.argv[1] == 'heatmap':
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --otu_table     Input OTU table (Required)
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
            sys.exit
    
    elif sys.argv[1] == 'taxonomy':
        #look in DB folder for databses
        db_list = ['DB_name', 'FASTA originated from', 'Fwd Primer', 'Rev Primer', 'Records']
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
        if len(db_list) < 6:
            db_print = "No DB configured, run 'ufits database' command for format database."
        else:
            d = flatten(db_list)
            db_print = fmtcols(d, 5)
        
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --fasta         Input FASTA file (i.e. OTUs from ufits cluster) (Required)
            -o, --out           Base name for output file. Default: ufits-taxonomy.<method>.txt
            -m, --method        Taxonomy method. Default: utax [utax, usearch] (Required)
            -d, --db            Database (must be in UDB format).
            --append_taxonomy   OTU table to append taxonomy. Default: none
            --utax_cutoff       UTAX confidence value cutoff. Default: 0.8 [0 to 0.9]
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
            try:
                dbLocation = arguments.index('-d')
            except ValueError:
                dbLocation = arguments.index('--db')
            dbLocation = dbLocation + 1
            if arguments[dbLocation] in okay_list:
                arguments[dbLocation] = os.path.join(script_path, 'DB', arguments[dbLocation])
                subprocess.call(arguments)
            else:
                print "Database %s not found, please choose a Configured DB" % arguments[dbLocation]
                print help
                sys.exit
        else:
            print help
            sys.exit
                    
    elif sys.argv[1] == 'database':
        help = """
Usage:      ufits %s <arguments>
version:    %s

Description:    Setup/Format reference database for ufits taxonomy command.
    
Arguments:  -i, --fasta         Input FASTA file (UNITE DB or UNITE+INSDC)
            -o, --out           Base Name for Output Files. Default: DB/ufits
            -f, --fwd_primer    Forward primer. Default: GTGARTCATCGAATCTTTG (fITS7)
            -r, --rev_primer    Reverse primer. Default: TCCTCCGCTTATTGATATGC (ITS4)
            --unite2utax        Reformat FASTA headers to UTAX format. Default: on
            --drop_ns           Removal sequences that have > x N's. Default: 8
            --create_db         Create a DB. Default: usearch [utax, usearch]
            --skip_trimming     Keep full length sequences. Default: off (not recommended)
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
            sys.exit
    
    elif sys.argv[1] == 'download':
        help = """
Usage:      ufits %s <arguments>
version:    %s

Description:    Download reference databases for ufits taxonomy command.
    
Arguments:  -i, --input     Database to download. [unite, unite_insd, rtl_ITS, all]
     
            
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
            sys.exit
    
    elif sys.argv[1] == 'summarize':
        help = """
Usage:      ufits %s <arguments>
version:    %s
    
Arguments:  -i, --table     OTU Table containing Taxonomy information (Required)
            -o, --out       Base name for output files. Default: ufits-summary
            --format       Image output format. Default: eps [eps, svg, png, pdf]
            --percent       Convert numbers to Percent of Sample. Default: off
            
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
            sys.exit
    elif sys.argv[1] == 'version' or '-version' or '--version':
        print "ufits v.%s" % version
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        sys.exit()
    
    
else:
    print default_help
        