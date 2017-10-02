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

# get arguments and load design
args <- commandArgs(trailingOnly=T)
design <- read.table("design_all_samples.txt", header=T, stringsAsFactors=F, sep="\t")
tissues <- design$tissue

# define vector of which samples are to be used,
# i.e. are matching the keywords
use <- rep(F, nrow(design))
ofile <- "sample_definitions/"
for(i in 1:length(args)){
  a <- args[i]
  a <- tolower(strsplit(a, " ")[[1]])
  u <- NULL
  for(p in a){
    if(is.null(u)){
      u <- grepl(p, tissues, ignore.case=T)
    } else {     
      u <- u & grepl(p, tissues, ignore.case=T)
    }

    # first addition, avoid the '_' 
    if(nchar(ofile)==19){
      ofile <- paste0(ofile, p)
    } else {
      ofile <- paste0(ofile, "_", p)
    }
  }
  use <- use | u
}

# add file extension
ofile <- paste0(ofile, ".R")

# get the matching samples and write to output file
# (as a R-vector definition)
samples <- design[use,]$sample
cat("Found", length(samples), "samples.\n")
cat(file=ofile, "samp <- c(")
for(i in 1:length(samples)) {
  if(i==length(samples)){
    cat(file=ofile, paste0("\"", samples[i], "\")"), append=T)
  } else {
    cat(file=ofile, paste0("\"", samples[i], "\",\n"), append=T)
  }
}
