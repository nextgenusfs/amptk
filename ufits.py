#!/usr/bin/env python

#Wrapper script for UFITS package.

import sys, os, subprocess, inspect
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

version = '0.1.1'

default_help = """
Usage:      ufits <command> <arguments>
version:    %s
    
Command:    ion         pre-process Ion Torrent data (find barcodes, trim/pad)
            illumina    pre-process folder of de-multiplexed Illumina data (gunzip, merge PE, trim/pad, concatenate)
            cluster     cluster OTUs (using UPARSE algorithm)
            filter      OTU table filtering
            heatmap     Create heatmap from OTU table
            
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
            cmd = os.path.join(script_path, 'ufits-process_ion.py')
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
            -f, --fwd_primer    Forward primer sequence. Default: GTGARTCATCGAATCTTTG (fITS7)
            -r, --rev_primer    Reverse primer sequence Default: TCCTCCGCTTATTGATATGC (ITS4)
            -n, --name_prefix   Prefix for re-naming reads. Default: R_
            -m, --min_len       Minimum length read to keep. Default: 50
            -l, --trim_len      Length to trim/pad reads. Default: 250
            -u, --usearch       USEARCH executable. Default: usearch8
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'ufits-process_illumina_folder.py')
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
            cmd = os.path.join(script_path, 'ufits-OTU_cluster.py')
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
            -p, --index_bleed   Filter index bleed between samples (percent). Default: 0.1
            -d, --delimiter     Delimiter of OTU table. Default: tsv  [csv, tsv] 
            -n, --names         Change names of barcodes (CSV mapping file). Default: off
            --mc                Mock community FASTA file. Default: ufits_mock3.fa
            --col_order         Column order (comma separated list). Default: sort naturally
            --convert_binary    Convert OTU table to binary (1's and 0's). Default: off
            --trim_data         Filter the data. Default: on [on, off]
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com      
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'ufits-mock_filter.py')
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
            -d, --delimiter     Delimiter of OTU table. (Required) [csv, tsv]    
            --format            Image output format. Default: eps [eps, svg, png, pdf]
            --col_order         Column order (comma separated list). Default: sort naturally
            --square            Maintain aspect ratio. Default: off
            --colors            Color Palette. Default: Greys [Blues, Reds, Purples, Greens, BuPu, BuGn, OrRd, YlGnBu, YlOrRd]
            --percent           Convert numbers to Percent of Sample. Default: off
            --min_num           Minimum reads per OTU to graph. Default: 1
            
Written by Jon Palmer (2015) nextgenusfs@gmail.com   
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'ufits-heatmap.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        sys.exit()
    
else:
    print default_help
        