#!/usr/bin/env Rscript

#run like Rscript --vanilla data2_pipeline.R input_folder output.csv platform pool cpus
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
install.packages.auto("dada2")
install.packages.auto("ShortRead")
#install.packages.auto("ggplot2")

#setup arguments, first is path to folder, second is path to output csv
args = commandArgs(trailingOnly=TRUE)

#load packages
library(dada2)
library(ShortRead)

#print the packages and versions
print("-------------")
print(paste("R", getRversion()))
for (package_name in sort(loadedNamespaces())) {
    print(paste(package_name, packageVersion(package_name)))
}
print("-------------")
print("Loading Data from folder")

#get DADA2 version
DADAversion <- packageVersion("dada2")
if (DADAversion >= '1.1.1') {
    CORES <- as.integer(args[5])
} else {
    CORES <- 'FALSE'
}

#load the data from a folder
path <- args[1]
fns <- list.files(path)
Seqs <- file.path(path, fns)
Seqs

#get sample names
sample.names <- sapply(strsplit(fns,".fastq"), `[`,1)
sample.names

#Dereplication
print("-------------")
print("Dereplicating data")
derepSeqs <- derepFastq(Seqs, verbose=TRUE)

#name the derep class with sample names
names(derepSeqs) <- sample.names

#Sample inference
print("-------------")
print("Sample inference")
if (args[3] == 'illumina') {
    dadaSeqs <- dada(derepSeqs, err=NULL, selfConsist=TRUE, pool=args[4], USE_QUALS=TRUE, multithread=CORES)
} else if (args[3] == 'ion') {
    dadaSeqs <- dada(derepSeqs, err=NULL, selfConsist=TRUE, pool=args[4], HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, USE_QUALS=TRUE, multithread=CORES)
} else {
    dadaSeqs <- dada(derepSeqs, err=NULL, selfConsist=TRUE, pool=args[4], HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, USE_QUALS=FALSE, multithread=CORES)
}

#make sequence table
seqtab <- makeSequenceTable(dadaSeqs, orderBy = "abundance")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

#transpose
transTable <- t(seqtab.nochim)

#output sequence table to csv
write.csv(transTable, file=args[2])
