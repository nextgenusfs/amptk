#!/usr/bin/env Rscript
#script to run LULU on OTU table in R

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

is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
  } 

#check for LULU package, install if not there
if (!is.installed('lulu')) {
    install.packages.auto("dplyr")
    install.packages.auto("devtools")
    install_github("tobiasgf/lulu")
}

#import necessary packages
library("lulu"); packageVersion("lulu")

#setup arguments, first is OTU table, second is match table, 3rd thru 6th are options, and 7th is output OTU table, 8th list of OTUs to keep, 9th is mapping file
args = commandArgs(trailingOnly=TRUE)
print(args)

#load in the OTU table and match file
otutab <- read.csv(args[1], sep='\t', header=TRUE, as.is=TRUE, row.names = 1)
matchlist <- read.table(args[2], header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

#now run the script
#curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
cr <- lulu(otutab, matchlist, minimum_ratio_type = args[3], minimum_ratio = args[4], minimum_match = args[5], minimum_relative_cooccurence = args[6])

#sort the resulting OTU table
t <-data.frame("OTUID"=rownames(cr$curated_table), cr$curated_table)
sortT <- t[with(t, order(as.integer(sub('\\D+', '', OTUID)))),]
colnames(sortT)[which(names(sortT) == "OTUID")] <- "#OTU ID"
write.table(sortT, file=args[7], sep='\t', row.names=FALSE, quote=FALSE)
write(cr$discarded_otus, file=args[8])
write.table(cr$otu_map, file=args[9], sep='\t', quote=FALSE)

