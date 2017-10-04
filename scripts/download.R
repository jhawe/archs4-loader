###################################################################
#' The code has been adapted from the code given at               #
#' http://amp.pharm.mssm.edu/archs4/help.html                     #
#'                                                                #
#' We essentially extract the needed data from the downloaded     #
#' archive and normalize it if necessary. Additionally we create  #
#' some basic diagnostic plots on the data.                       #
#'                                                                #
#' @author Johann Hawe                                            #
#'                                                                #
###################################################################

# R script to download selected samples
# Copy code and run on a local machine to initiate download

# Check for dependencies and install if missing
packages <- c("rhdf5", "preprocessCore", "sva")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    print("Install required packages")
    source("https://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
    biocLite("preprocessCore")
    biocLite("sva")
}
library("rhdf5")
library("preprocessCore")
library("sva")

# check package for option parsing
if(!require(optparse)){
  install.packages(optparse)
}

# parse arguments
# get arguments and define global vars
opt = parse_args(OptionParser(option_list=list(
  make_option(c("-s", "--samples"), type="character", default=NULL),
  make_option("--normalize", type="logical", action="store_true", default=FALSE),
  make_option("--sva", type="logical", action="store_true", default=FALSE),
  make_option("--peer", type="logical", action="store_true", default=FALSE))));

# define global vars
samp.file = opt$samples
if(!file.exists(samp.file)){
  stop("Specified sample-file does not exist!")
}
NORM = opt$normalize
SVA = opt$sva
PEER = opt$peer
if(SVA & PEER) {
  stop("Choose either COMBAT or PEER!")
}
if(PEER) {
  source("scripts/peer.R")
}

# get the samples
cat("Using samples defined in:", samp.file, "\n")
source(samp.file)
keyword <- gsub("\\..*", "", basename(samp.file))
dir.create(keyword)
print(paste0("Saving data under ./", keyword, "/"))

# the actual archive file with data
destination_file = "human_matrix_download.h5"

# name of output expression file
# depending in processing type
`%+%` = paste0
f = "expression_matrix"
if(NORM) {
  f = f %+% "_norm.tsv"
}
# peer or sva normalization?
if(PEER) {
  f = f %+% "_peer.tsv"
} else if(SVA) {
  f = f %+% "_sva.tsv"
}
extracted_expression_file = paste0(keyword, "/", f)

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
# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# get info for selected samples
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
des = des[sample_locations,,drop=F]

# write design file
write.table(des, file=file.path(keyword, "design.txt"), sep="\t", quote=T, col.names=T, row.names=F)
print(paste0("Design file was created at ", getwd(), "/", keyword, "/design.txt"))
rm(des)

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()

# normalize samples and correct for differences in gene count distribution
if(NORM | SVA | PEER){
  print("Normalizing data.")
  expression = log2(expression+1)
  expression = normalize.quantiles(expression)
}

rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# correct batch effects in gene expression
if(SVA) {
  print("Removing batch effects using ComBat.")
  series = series[sample_locations]
  batchid = match(series, unique(series))
  expression <- ComBat(dat=expression, batch=batchid, par.prior=TRUE, prior.plots=FALSE)
} else if(PEER) {
  print("Removing batch effects using PEER.")
  # peer expects an NxG matrix (N=#samples)
  expression <- t(correct.peer(data=t(expression)))
}

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE, col.names=NA, row.names=T)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

# Create histogram of expression values and heatmap
# of gene correlations
print("Creating heatmap and expression histogram.")

pdf(gsub("\\.tsv$", ".pdf", extracted_expression_file))

# boxplot and histogram of max 300 random samples
toplot <- expression[,sample(1:ncol(expression),min(ncol(expression), 300))]
boxplot(as.data.frame(toplot), main="expression of random samples", xlab="samples", ylab="expression")
hist(toplot, breaks=100, main="expression values",xlab="expression")

# gene against gene correlation plots
cors <- cor(t(toplot))
hist(cors, breaks=100, main="correlation of genes", xlab="correlation")
dev.off()
