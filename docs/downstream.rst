
.. _downstream:

AMPtk Downstream Tools
================

AMPtk to Phyloseq
-------------------------------------
Importing a BIOM table from AMPtk into Phyloseq is relatively straightforward:

.. code-block:: R

    #load in biom and tree file
    physeq <- import_biom('mydata.biom', 'mydata.tree.phy')
    
    #rename taxonomy headings
    colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", 
      "Order", "Family", "Genus", "Species")

    #print the sample variables
    sample_variables(physeq)

    

AMPtk to QIIME
-------------------------------------
The BIOM file produced by ``amptk taxonomy`` is compatible directly with QIIME community ecology scripts, assuming that you have a QIIME compatible mapping file:

.. code-block:: none

    summarize_taxa_through_plots.py -i otu_table.biom -o taxa_summary -m mapping_file.txt
    

AMPtk to Online Tools
-------------------------------------
The BIOM output is also compatible with `Phinch <http://phinch.org>`_ web based visualization as well as `MetaCoMET <https://probes.pw.usda.gov/MetaCoMET/>`_.

