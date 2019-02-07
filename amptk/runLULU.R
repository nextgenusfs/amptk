#!/usr/bin/env Rscript
#script to run LULU on OTU table in R

#import necessary packages
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("lulu"); packageVersion("lulu")

#setup arguments, first is OTU table, second is match table, 3rd thru 6th are options, and 7th is output OTU table, 8th list of OTUs to keep, 9th is mapping file
args = commandArgs(trailingOnly=TRUE)
print(args)

#load in the OTU table and match file
otutab <- read.csv(args[1], sep='\t', header=TRUE, as.is=TRUE, row.names = 1, check.names=FALSE)
orig.cols <- colnames(otutab)
output.cols <- append(c('#OTU ID'), orig.cols)
colnames(otutab) <- make.names(colnames(otutab))
matchlist <- read.table(args[2], header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

#now run the script
#curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
cr <- lulu(otutab, matchlist, minimum_ratio_type = args[3], minimum_ratio = args[4], minimum_match = args[5], minimum_relative_cooccurence = args[6])

#sort the resulting OTU table
t <-data.frame("OTUID"=rownames(cr$curated_table), cr$curated_table)
sortT <- t[with(t, order(as.integer(sub('\\D+', '', OTUID)))),]
write.table(sortT, file=args[7], sep='\t', row.names=FALSE, quote=FALSE, col.names=output.cols)
write(cr$discarded_otus, file=args[8])
write.table(cr$otu_map, file=args[9], sep='\t', quote=FALSE)

