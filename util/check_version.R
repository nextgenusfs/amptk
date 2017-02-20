#!/usr/bin/env Rscript
Rversion <- R.Version()$version.string
dadaversion <- packageVersion("dada2")
parts <- strsplit(Rversion, ' ')
Rvers <- parts[[1]][3]
output <- paste(Rvers,dadaversion, sep=',')
cat(output,"\n")