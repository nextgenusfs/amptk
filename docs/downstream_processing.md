###UFITS to QIIME
Running the ["core diversity"](http://qiime.org/scripts/core_diversity_analyses.html) script in QIIME can be accomplished like so from UFITS output.
```
#run biom summary to see find lowest number of reads
biom summarize-table -i ufits.output.biom

#run core diversity script from QIIME
core_diversity_analyses.py -i ufits.output.biom -m ufits.mapping_file.txt -t ufits.tree.phy \
                           -c Treatment1,Treatment2 -e 90000 -o core_diversity_output
```

###UFITS to PhyloSeq
It is relatively straightforward to import your data into PhyloSeq in R.
```
#load packages
library("phyloseq")
library("ggplot2")
library("plyr")

#import biom and tree file into PhyloSeq
physeq <- import_biom('/path/to/ufits.output.biom', treefilename='/path/to/ufits.tree.phy')

#for example you can now run alpha diversity on your samples
plot_richness(physeq, measures=c("Observed", "Chao1", "Shannon"))

#to group your samples based on treatment
plot_richness(physeq, x=Treatment1, color=Treatment1, measures=c("Observed", "Chao1", "Shannon"))
```
