#!/bin/usr/Rscript

###################################################################
#'                                                                #
#' The scripts parses the available design file for all available #
#' samples and greps in the tissue column for the specified       #
#' keywords. Keywords enclosed by quotation marks are connected   #
#' with 'AND', otherwise 'OR' connection is performed.            #
#'                                                                #
#' @author Johann Hawe                                            #
#'                                                                #
###################################################################

# load libraries
library(tidyverse)

# get arguments and load design
design <- read_tsv(snakemake@input$design)
tissues <- design$tissue

# get output file
fout <- snakemake@output[[1]]

# get keywords
keywords <- snakemake@wildcards$keywords
# we assume "_" in the keywords for word separation
keywords <- strsplit(keywords, "_")[[1]]

# define vector of which samples are to be used,
# i.e. are matching the keywords
use <- NULL
for(k in keywords){
  if(is.null(use)){
    use <- grepl(k, tissues, ignore.case=T)
  } else {     
    use <- use & grepl(k, tissues, ignore.case=T)
  }
}

# get the matching samples and write to output file
# (as a R-vector definition)
samples <- design[use,]$sample
print(paste0("Found ", length(samples), " samples."))
cat(file=fout, "samp <- c(")
for(i in 1:length(samples)) {
  if(i==length(samples)){
    cat(file=fout, paste0("\"", samples[i], "\")"), append=T)
  } else {
    cat(file=fout, paste0("\"", samples[i], "\",\n"), append=T)
  }
}
