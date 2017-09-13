
.. _quick-start:

AMPtk Quick Start
================================

AMPtk can process NGS amplicon data from either 454, Illumina, or Ion Torrent sequencing platforms.
In order to get started, you will need to know how your amplicon reads were generated, i.e.
which primers you have used, PE or SE reads, where the amplicons have been barcoded, etc. For more
information see the :ref:`ReadLayout` information.

All AMPtk commands show a help menu when they are run without any parameters.

AMPtk comes with some very small test datasets, here we will go through a standard MiSeq Illumina test run. The data is located in the ``test_data`` folder of the amptk installation directory.  Move to that folder to begin.

Pre-Processing Reads
-------------------------------------
The first step is to pre-process or demultiplex your reads.  In this example, the data is PE MiSeq (2 x 300 bp) and it has been processed by the sequencing center using the Illumina program ``bcl2fastq`` so each sample has a forward and reverse read file (denoted by _R1 and _R2 in the sample name). The folder ``illumina_test_data`` contains 3 PE samples we will use for a demo.  Since the data is in this format, we will first use the ``amptk illumina`` format.  Note that some sequencing centers distribute their data in different formats, ``amptk illumina2`` and ``amptk illumina3`` are additional methods that can be used to pre-process illumina data.

.. code-block:: none

    #run the pre-processing script, specifying primers used
    amptk illumina -i illumina_test_data -o miseq -f fITS7 -r ITS4

You should see that script found 600 valid reads, and output them into a file called ``miseq.demux.fq.gz`` and you'll note that it also output a QIIME-like mapping file ``miseq.mapping_file.txt`` which can be used later on to add meta data to and generate a BIOM file (during assigning taxonomy step).  You can find a detailed log file in ``miseq.amptk-demux.log``.  The script has merged the PE reads using usearch, removed phiX spike-in, removed the forward and reverse primers, relabeled the FASTQ headers, and then concatenated all samples together.  Note that the default settings require that reads have a valid forward primer, if you used a custom sequencing primer then you will need to turn this setting off with ``--require_primer off``.

Clustering Data
-------------------------------------

