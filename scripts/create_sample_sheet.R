#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#' The code has been adapted from the code given at
#' http://amp.pharm.mssm.edu/archs4/help.html
#'
#' We simply use the main h5 file to get basic information for
#' all the samples.
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------

library("rhdf5")

source("scripts/lib.R")

# name of output design file
design_file = snakemake@output[[1]]

# the actual archive file containing the data
destination_file = snakemake@output[[2]]

# Check if gene expression file was already downloaded, if it doesn't exist
# then download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
    download.file(url, destination_file, quiet = FALSE)
} else{
    print("Local data file already exists.")
}

design <- load_design(destination_file)

# write design file
write.table(design, file=design_file, sep="\t",
  quote=T, col.names=T, row.names=F)

# done
sessionInfo()
