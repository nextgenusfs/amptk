###AMPtk to QIIME
Running the ["core diversity"](http://qiime.org/scripts/core_diversity_analyses.html) script in QIIME can be accomplished like so from AMPtk output.
```
#run biom summary to see find lowest number of reads
biom summarize-table -i amptk.output.biom

#run core diversity script from QIIME
core_diversity_analyses.py -i amptk.output.biom -m amptk.mapping_file.txt -t amptk.tree.phy \
                           -c Treatment1,Treatment2 -e 90000 -o core_diversity_output
```

###AMPtk to PhyloSeq
It is relatively straightforward to import your data into PhyloSeq in R.
```
#load packages
library("phyloseq")
library("ggplot2")
library("plyr")

#import biom and tree file into PhyloSeq
physeq <- import_biom('/path/to/amptk.output.biom', treefilename='/path/to/amptk.tree.phy')

#for example you can now run alpha diversity on your samples
plot_richness(physeq, measures=c("Observed", "Chao1", "Shannon"))

#to group your samples based on treatment
plot_richness(physeq, x=Treatment1, color=Treatment1, measures=c("Observed", "Chao1", "Shannon"))
```

###AMPtk to Vegan
Here you can take adavantage of the biom package to load your data and easily convert it to a format that vegan requires (which is essentially a classic OTU table that has been transposed).
```
#load packages
library("biom")
library("vegan")

#load biom file c
b1 <- read_biom('/path/to/amptk.output.biom')

#split metadata, taxonomy, and otu_table
otu_table <- t(as.data.frame(as.matrix(biom_data(b1))))
metadata <- sample_metadata(b1)
taxonomy <- observation_metadata(b1)

#convert otu table to binary
otu_table[otu_table>0] <-1

#create distance matrix for binary data
rp <- raupcrick(otu_table, null="r1", nsimul=999, chase=FALSE)

#run ordination
metaMDS(rp, k=2, autotransform=FALSE, wascores=F, trymax=200)

#run hypothesis adonis test on Treatment1 variable
adonis(rp~metadata$Treatment1, permutations=9999)
```