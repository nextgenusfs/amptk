#!/usr/bin/env Rscript

country.code <- 'us'  # use yours
url.pattern <- 'https://'  # use http if you want
repo.data.frame <- subset(getCRANmirrors(), CountryCode == country.code & grepl(url.pattern, URL))
options(repos = repo.data.frame$URL)

#check for required packages, install if necessary
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    #update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

install.packages.auto("phyloseq")
install.packages.auto("ggplot2")
install.packages.auto("plyr")
install.packages.auto("vegan")
install.packages.auto("RColorBrewer")
install.packages.auto("indicspecies")

#experimenting with phyloseq, import necessary packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("vegan"); packageVersion("vegan")
library("RColorBrewer"); packageVersion("RColorBrewer")
library("indicspecies"); packageVersion("indicspecies")

#setup arguments, first is path to folder, second is path to output csv
args = commandArgs(trailingOnly=TRUE)
print(args)
lenArgs <- length(args)
if ( lenArgs > 5) {
    dropOTUs <- args[6:lenArgs]
} 

#functions
#convert phyloseq object to vegan matrix
veganotu = function(physeq) {
    require("vegan")
    OTU = otu_table(physeq)
    if (taxa_are_rows(OTU)) {
        OTU = t(OTU)
    }
    return(as(OTU, "matrix"))
}

#function to draw phyloseq object into NMDS ordination with centroids and standard deviation error bars
nmds_pretty <- function(physeq, ord, dataset, variable, colors, heading, custom_theme) { 
    scrs <- scores(ord)
    centroids <- setNames(aggregate(scrs, by=list(dataset[[variable]]), FUN=mean), c(variable, "NMDS1", "NMDS2"))
    stdev <- setNames(aggregate(scrs, by=list(dataset[[variable]]), FUN=sd), c(variable, "sd1", "sd2"))
    combined <- setNames(merge(centroids, stdev, by=variable), c(variable, "NMDS1", "NMDS2", "sd1", "sd2"))
    p <- plot_ordination(physeq, ord, type="samples", color=variable)
    p <- p + ggtitle(paste(heading)) + theme(plot.title = element_text(hjust = 0.5)) + custom_theme
    p <- p + xlim(min(scrs[,1]),max(scrs[,1])) + ylim(min(scrs[,2]),max(scrs[,2]))
    p <- p + geom_point(data=combined, size=4)
    p <- p + geom_errorbar(data=combined, aes(ymin=NMDS2-sd2, ymax=NMDS2+sd2), width=0.1) 
    p <- p + geom_errorbarh(data=combined, aes(xmin=NMDS1-sd1, xmax=NMDS1+sd1), height=0.1)
    p <- p + scale_color_manual(values=colors)
    return(p)
}

#load in biom and tree file
physeqLoad <- import_biom(args[1], treefilename=args[2])

#get number of taxonomic ranks from phyloseq object
numRanks <- length(rank_names(physeqLoad))
if ( numRanks > 7 ) {
	colnames(tax_table(physeqLoad)) <- c("Domain", "Kingdom", "Phylum", "Class", 
  "Order", "Family", "Genus", "Species")
  } else {
	colnames(tax_table(physeqLoad)) <- c("Kingdom", "Phylum", "Class", 
  "Order", "Family", "Genus", "Species")  
}

#check if OTUs to drop
if ( lenArgs > 5 ) {
    goodTaxa <- setdiff(taxa_names(physeqLoad), dropOTUs)
    physeq <- prune_taxa(goodTaxa, physeqLoad)
} else {
    physeq <- physeqLoad 
}

#print the sample variables
sample_variables(physeq)

#get vegan OTU table
votu <- veganotu(physeq)

#make sure vegan OTU is binary (this is necessary for RaupCrick)
votu[votu>0] <-1

#recreate binary table as phyloseq object
OTU = otu_table(t(votu), taxa_are_rows = TRUE)
physeq2 = phyloseq(OTU, tax_table(physeq), sample_data(physeq), phy_tree(physeq))

#loop through treatments, remove the non treatment options in biom that are relic of 
remove <- c("LinkerPrimerSequence","phinchID","RevBarcodeSequence","DemuxReads","BarcodeSequence","ReversePrimer")
treatments <- setdiff(sample_variables(physeq), remove)
print(treatments)

#all white theme
t1 <-theme(                              
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  panel.border = element_rect(color = "black", fill=NA, size=0.4)
)

binary_distances <- c("unifrac", "jaccard")
special_distances <- c("raupcrick")
abundance_distances <- c("bray", "wunifrac")


#loop through all each of the treatments
for (y in 1:length(treatments)) {
    #first need to check if for this treatment we have no_data, if so need to drop
    print(treatments[y])
    variables = get_variable(physeq, treatments[y])
    num_variables = length(unique(variables))
    if ( is.element('no_data',variables) ) {
        num_variables <- num_variables - 1
        }
    if ( num_variables > 1 ) {
        remove_idx = as.character(get_variable(physeq, treatments[y])) != "no_data"
        PhyBray = prune_samples(remove_idx, physeq)
        PhyRaup = prune_samples(remove_idx, physeq2)
        print(PhyBray)
        
        #get metadata again for updated phyloseq object
        metadata <- as(sample_data(PhyBray), "data.frame")
    
        #convert to vegan OTU matrix for hypothesis test
        votuRaup <- veganotu(PhyRaup)
    
        #run distance for raupcrick
        if ( is.element(args[4], special_distances) ) {
            dist <- raupcrick(votuRaup, null="r1", nsimul=999, chase=FALSE)
        } else if ( is.element(args[4], binary_distances) ) {
            dist <- distance(PhyRaup, method=args[4])
        } else if ( is.element(args[4], abundance_distances) )  {
            dist <- distance(PhyBray, method=args[4])
        } else 
            print("Distance not supported")
        
        #setup output stats file
        statsout <- paste(paste(dirname(args[1]), args[3], args[3], sep='/'),'.', treatments[y],'.adonis-', args[4], '.txt', sep='')
    
        #hypothesis test for significance
        testA <- adonis(dist~metadata[[treatments[y]]], permutations=9999)
        betaA <- betadisper(dist, metadata[[treatments[y]]])
        pA <- permutest(betaA, permutations=9999)
        capture.output(testA,file=statsout)
        capture.output(pA,file=statsout, append=TRUE)

        #setup color palette, if more than 15 don't even draw it as it is dumb
        num_colors = length(unique(get_variable(PhyBray, treatments[y])))
        print(num_colors)
        getPalette = colorRampPalette(brewer.pal(9, "Set1"))
        if (num_colors < 10 ) {
            colors = brewer.pal(9, "Set1")
        } else if ( num_colors > 15) {
            next
        } else 
            colors = getPalette(num_colors)
        
        #run ordination on previously computed distances, then save in a 2x2 image with shared legend.  
        ord <- ordinate(PhyRaup, method=args[5], distance=dist)
        capture.output(ord, file=statsout, append=TRUE)
        p <- nmds_pretty(PhyRaup, ord, metadata, treatments[y], colors, args[4], t1)
        ggsave(paste(paste(dirname(args[1]),args[3], args[3], sep='/'),'.', treatments[y],'.',args[4],'.ordination.pdf', sep=''), plot=p)
    
        #alpha diversity
        pr <- plot_richness(PhyBray, x=treatments[y], color=treatments[y], measures=c("Observed", "Chao1", "Shannon"))
        pr = pr + geom_point(size=3, alpha=0.4)
        pr = pr + scale_color_manual(values=colors)
        ggsave(paste(paste(dirname(args[1]),args[3], args[3], sep='/'),'.', treatments[y],'.alpha_diversity.pdf', sep=''), plot=pr)
    
        #run indicator species test
        indicsummary <- paste(paste(dirname(args[1]), args[3], args[3], sep='/'),'.', treatments[y],'.indicator_species_summary.txt', sep='')
        indicresults <- paste(paste(dirname(args[1]), args[3], args[3], sep='/'),'.', treatments[y],'.indicator_species_results.txt', sep='')
        if ( is.element(args[4], abundance_distances) ) {
            indspVOTU <- as.data.frame(veganotu(PhyBray))
        } else 
            indspVOTU <- as.data.frame(votuRaup)
        indicval <- multipatt(indspVOTU, metadata[[treatments[y]]], control=how(nperm=999))
        capture.output(summary(indicval, indvalcomp = TRUE), file=indicsummary)
        capture.output(coverage(indspVOTU, indicval), file=indicsummary, append=TRUE)
        indicAll <- indicval$sign
        indicAll.pvalues <- indicAll[,"p.value"]
        pvalue.fdr.correct <- p.adjust(indicAll.pvalues, method="fdr")
        indicAll.fdr <- cbind(indicAll, pvalue.fdr.correct)
        attach(indicAll.fdr)
        indic.fdr.sort <- indicAll.fdr[order(pvalue.fdr.correct, p.value, na.last=NA),]
        capture.output(indic.fdr.sort, file=indicresults)
        detach(indicAll.fdr)
        #now run indicator species for correlation indices, i.e. phi from the docs
        corr_summary <- paste(paste(dirname(args[1]), args[3], args[3], sep='/'),'.', treatments[y],'.indicator_correlation_summary.txt', sep='')
        corr_results <- paste(paste(dirname(args[1]), args[3], args[3], sep='/'),'.', treatments[y],'.indicator_correlation_results.txt', sep='')
        phi <- multipatt(indspVOTU, metadata[[treatments[y]]], func="r.g", control=how(nperm=999))
        capture.output(summary(phi, indvalcomp = TRUE), file=corr_summary)
        capture.output(coverage(indspVOTU, phi), file=corr_summary, append=TRUE)
        phiAll <- phi$sign
        phiAll.pvalues <- phiAll[,"p.value"]
        pvalue.fdr.correct <- p.adjust(phiAll.pvalues, method="fdr")
        phiAll.fdr <- cbind(phiAll, pvalue.fdr.correct)
        attach(phiAll.fdr)
        phi.fdr.sort <- phiAll.fdr[order(pvalue.fdr.correct, p.value, na.last=NA),]
        capture.output(phi.fdr.sort, file=corr_results)
        detach(phiAll.fdr)       
        
    } else
        print("skipping")
        next             
}
