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
install.packages.auto("plotly")
install.packages.auto("htmltools")
install.packages.auto("DT")


#experimenting with phyloseq, import necessary packages
library("plyr"); packageVersion("plyr")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan"); packageVersion("vegan")
library("RColorBrewer"); packageVersion("RColorBrewer")
library("plotly"); packageVersion("plotly")
library("htmltools"); packageVersion("htmltools")
library("DT"); packageVersion("DT")

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

#richness filter to drop samples if less than n taxa
phyloseq_richness_filter <- function(physeq, mintaxa = 10){
  sp <- estimate_richness(physeq, measures = "Observed")
  samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]
  if(length(samples_to_keep) == 0){
    stop("All samples will be removed.\n")  
  }
  if(length(samples_to_keep) == nsamples(physeq)){
    cat("All samples will be preserved\n")
    res <- physeq
  }
  if(length(samples_to_keep) < nsamples(physeq)){
    removed <- setdiff(sample_names(physeq), samples_to_keep)
    cat("Removing", length(removed), "samples with fewer than", mintaxa, "taxa:", removed, "\n")
    res <- filter_taxa(prune_samples(samples = samples_to_keep, x = physeq), function(x) sum(x) > 1, TRUE)
  }
  return(res)
}

#function to draw phyloseq object into NMDS ordination with centroids and standard deviation error bars
nmds_pretty <- function(physeq, ord, dataset, variable, colors, heading, custom_theme, centroids = TRUE) { 
    p <- plot_ordination(physeq, ord, type="samples", color=variable)
    p <- p + geom_point(aes(text=sprintf("ID: %s", sample_names(physeq))))
    p <- p + ggtitle(paste(heading)) + theme(plot.title = element_text(hjust = 0.5)) + custom_theme
    if (centroids == TRUE) {
      scrs <- scores(ord)
      centroids <- setNames(aggregate(scrs, by=list(dataset[[variable]]), FUN=mean), c(variable, "NMDS1", "NMDS2"))
      stdev <- setNames(aggregate(scrs, by=list(dataset[[variable]]), FUN=sd), c(variable, "sd1", "sd2"))
      combined <- setNames(merge(centroids, stdev, by=variable), c(variable, "NMDS1", "NMDS2", "sd1", "sd2"))
      p <- p + xlim(min(scrs[,1]),max(scrs[,1])) + ylim(min(scrs[,2]),max(scrs[,2]))
      p <- p + geom_point(data=combined, size=4)
      p <- p + geom_errorbar(data=combined, aes(ymin=NMDS2-sd2, ymax=NMDS2+sd2), width=0.1) 
      p <- p + geom_errorbarh(data=combined, aes(xmin=NMDS1-sd1, xmax=NMDS1+sd1), height=0.1)
    } else {
      p <- p + stat_ellipse(type = "norm", linetype = 2)
      p <- p + stat_ellipse(type = "t")
    }
    p <- p + scale_color_manual(values=colors)
    return(p)
}

cca_pretty <- function(physeq, ord, dataset, variable, colors, heading, custom_theme) { 
    scrs <- scores(ord)
    centroids <- setNames(aggregate(scrs, by=list(dataset[[variable]]), FUN=mean), c(variable, "CCA1", "CCA2"))
    stdev <- setNames(aggregate(scrs, by=list(dataset[[variable]]), FUN=sd), c(variable, "sd1", "sd2"))
    combined <- setNames(merge(centroids, stdev, by=variable), c(variable, "CCA1", "CCA2", "sd1", "sd2"))
    p <- plot_ordination(physeq, ord, type="samples", color=variable)
    p <- p + ggtitle(paste(heading)) + theme(plot.title = element_text(hjust = 0.5)) + custom_theme
    p <- p + xlim(min(scrs[,1]),max(scrs[,1])) + ylim(min(scrs[,2]),max(scrs[,2]))
    p <- p + geom_point(data=combined, size=4)
    p <- p + geom_errorbar(data=combined, aes(ymin=NMDS2-sd2, ymax=NMDS2+sd2), width=0.1) 
    p <- p + geom_errorbarh(data=combined, aes(xmin=NMDS1-sd1, xmax=NMDS1+sd1), height=0.1)
    p <- p + scale_color_manual(values=colors)
    return(p)
}

#all white theme
t1 <-theme(                              
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  panel.border = element_rect(color = "black", fill=NA, size=0.4)
)

#function to output html page
save_tags <- function (tags, file, selfcontained = F, libdir = "./lib") 
{
  if (is.null(libdir)) {
    libdir <- paste(tools::file_path_sans_ext(basename(file)), 
                    "_files", sep = "")
  }
  htmltools::save_html(tags, file = file, libdir = libdir)
  if (selfcontained) {
    if (!htmlwidgets:::pandoc_available()) {
      stop("Saving a widget with selfcontained = TRUE requires pandoc. For details see:\n", 
           "https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md")
    }
    htmlwidgets:::pandoc_self_contained_html(file, file)
    unlink(libdir, recursive = TRUE)
  }
  return(htmltools::tags$iframe(src= file, height = "400px", width = "100%", style="border:0;"))
}


#######Start here#########

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

#remove samples if less than 3 taxa -- these likely to be outliers
physeqfilt <- phyloseq_richness_filter(physeq, mintaxa = 2)

#print the sample variables
sample_variables(physeqfilt)

#get vegan OTU table
votu <- veganotu(physeqfilt)

#make sure vegan OTU is binary (this is necessary for RaupCrick)
votu[votu>0] <-1

#recreate binary table as phyloseq object
OTU = otu_table(t(votu), taxa_are_rows = TRUE)
physeq2 = phyloseq(OTU, tax_table(physeqfilt), sample_data(physeqfilt), phy_tree(physeqfilt))

#loop through treatments, remove the non treatment options in biom that are relic of 
remove <- c("LinkerPrimerSequence","phinchID","RevBarcodeSequence","DemuxReads","BarcodeSequence","ReversePrimer")
treatments <- setdiff(sample_variables(physeqfilt), remove)
print(treatments)

binary_distances <- c("unifrac", "jaccard")
special_distances <- c("raupcrick", "aitchison")
abundance_distances <- c("bray", "wunifrac")


#loop through all each of the treatments
for (y in 1:length(treatments)) {
    #first need to check if for this treatment we have no_data, if so need to drop
    print(treatments[y])
    variables = get_variable(physeqfilt, treatments[y])
    num_variables = length(unique(variables))
    if ( is.element('no_data',variables) ) {
        num_variables <- num_variables - 1
        }
    if ( num_variables > 1 ) {
        remove_idx = as.character(get_variable(physeqfilt, treatments[y])) != "no_data"
        PhyBray <- filter_taxa(prune_samples(remove_idx, physeqfilt), function(x) sum(x) > 1, TRUE)
        PhyRaup <- filter_taxa(prune_samples(remove_idx, physeq2), function(x) sum(x) > 1, TRUE)
        print(PhyBray)
        
        #get metadata again for updated phyloseq object
        metadata <- as(sample_data(PhyBray), "data.frame")
    
        #convert to vegan OTU matrix for hypothesis test
        votuRaup <- veganotu(PhyRaup)
        votuBray <- veganotu(PhyBray)
    
        #run distance for raupcrick
        if ( is.element(args[4], special_distances) ) {
        	if ( args[4]=='raupcrick'){
        		d <- raupcrick(votuRaup, null="r1", nsimul=999, chase=FALSE)
        		PhyORD = PhyRaup
        	} else
        		d <- aDist(votuBray)
        		PhyORD = PhyBray
        } else if ( is.element(args[4], binary_distances) ) {
            d <- distance(PhyRaup, method=args[4])
            PhyORD = PhyRaup
        } else if ( is.element(args[4], abundance_distances) )  {
            d <- distance(PhyBray, method=args[4])
            PhyORD = PhyBray
        } else 
            print("Distance not supported")
            
        #hypothesis test for significance
        testA <- adonis(d~metadata[[treatments[y]]], permutations=9999)
        betaA <- betadisper(d, metadata[[treatments[y]]])
        pA <- permutest(betaA, permutations=9999)
        
        #create dataframe of stats from tests
        Df <- c(testA$aov.tab$Df[1], pA$tab$Df[1])
		Fstat <- c(testA$aov.tab$F.Model[1], pA$tab$F[1])
		R2 <- c(testA$aov.tab$R2[1], 'NA')
		Pvalue <- c(testA$aov.tab$`Pr(>F)`[1], pA$tab$`Pr(>F)`[1])
		xTab <- data.frame(Df,Fstat,R2,Pvalue, row.names=c('Adonis', 'Betadisper'))
		
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
        ord <- ordinate(PhyORD, method=args[5], distance=d)
        if ( args[5]== 'NMDS'){
        	p <- nmds_pretty(PhyORD, ord, metadata, treatments[y], colors, args[4], t1)
        } else
        	p <- cca_pretty(PhyORD, ord, metadata, treatments[y], colors, args[4], t1)
    
        #alpha diversity
        pr <- plot_richness(PhyBray, x=treatments[y], color=treatments[y], measures="Observed")
        pr <- pr + geom_point(size=3, alpha=0.4)
        pr <- pr + scale_color_manual(values=colors)
        
        #now create html of plotly outputs
        combined <- tags$div(
        	style = "display: flex; flex-wrap: wrap",
  			tags$div(ggplotly(p), style = "width: 90%; padding: 1em;"),
  			tags$div(ggplotly(pr), style = "width: 75%; padding: 1em;"),
  			tags$div(datatable(xTab), style = "width: 75%; padding: 1em;"))
        save_tags(combined, paste(paste(dirname(args[1]),args[3], args[3], sep='/'),'.', treatments[y],'.',args[4],'.html', sep=''))
        
    } else
        print("skipping")
        next             
}
