#!/usr/bin/env Rscript
is.installed <- function(mypkg){
is.element(mypkg, installed.packages()[,1])
} 

Rversion <- R.Version()$version.string
if (!is.installed("dada2")){
    dadaversion <- '0.0.0'
} else {
    dadaversion <- packageVersion("dada2")
    }
if (!is.installed("phyloseq")){
    phyloseqversion <- '0.0.0'
} else {
    phyloseqversion <- packageVersion("phyloseq")
    }
if (!is.installed("lulu")){
    luluversion <- '0.0.0'
} else {
    luluversion <- packageVersion("lulu")
    }
parts <- strsplit(Rversion, ' ')
Rvers <- parts[[1]][3]
output <- paste(Rvers, dadaversion, phyloseqversion, luluversion, sep=',')
cat(output,"\n")
