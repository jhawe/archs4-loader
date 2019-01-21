#!/usr/bin/Rscript

###################################################################
#' The code has been adapted from the code given at               #
#' http://amp.pharm.mssm.edu/archs4/help.html                     #
#'                                                                #
#' We simply use the main h5 file to get information for          #
#' all the samples.                                               #
#'                                                                #
#' @author Johann Hawe                                            #
#'                                                                #
###################################################################

# Check for dependencies and install if missing
packages <- c("rhdf5", "preprocessCore", "sva")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    print("Install required packages")
    source("https://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
}
library("rhdf5")

# the actual archive file with data
destination_file = "human_matrix_download.h5"

# name of output expression file
design_file = "design_all_samples.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
    download.file(url, destination_file, quiet = FALSE)
} else{
    print("Local data file already exists.")
}

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")
series = h5read(destination_file, "meta/Sample_series_id")
organism = h5read(destination_file, "meta/Sample_organism_ch1")
molecule = h5read(destination_file, "meta/Sample_molecule_ch1")
characteristics = h5read(destination_file, "meta/Sample_characteristics_ch1")
description = h5read(destination_file, "meta/Sample_description")
instrument = h5read(destination_file, "meta/Sample_instrument_model")

des = cbind(sample=samples, tissue, series, organism, molecule, characteristics, description,
               instrument)

# write design file
write.table(des, file=design_file, sep="\t", quote=T, col.names=T, row.names=F)

# done
sessionInfo()
