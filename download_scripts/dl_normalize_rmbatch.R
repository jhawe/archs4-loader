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

args <- commandArgs(trailingOnly=T)
if(length(args)==0){
  stop("You need to provide a sample file!")
}

# get the samples
samp.file <- args[1]
cat("Using samples defined in:", samp.file, "\n")
source(samp.file)
keyword <- gsub("\\..*", "", basename(samp.file))
dir.create(keyword)

destination_file = "human_matrix_download.h5"
extracted_expression_file = paste0(keyword, "/expression_matrix.tsv")

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
write.table(des, file=file.path(keyword, "design.txt"), sep="\t", quote=FALSE, col.names=T, row.names=F)
print(paste0("Design file was created at ", getwd(), "/", keyword, "/design.txt"))
rm(des)

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()

# normalize samples and correct for differences in gene count distribution
expression = log2(expression+1)
expression = normalize.quantiles(expression)

rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# correct batch effects in gene expression
series = series[sample_locations]
batchid = match(series, unique(series))
correctedExpression <- ComBat(dat=expression, batch=batchid, par.prior=TRUE, prior.plots=FALSE)

# Print file
write.table(correctedExpression, file=extracted_expression_file, sep="\t", quote=FALSE, col.names=NA, row.names=T)
print(paste0("Corrected expression file was created at ", getwd(), "/", extracted_expression_file))

# Create histogram of expression values and heatmap
# of gene correlations
print("Creating heatmap and expression histogram.")
library(pheatmap)
pdf(file.path(keyword, "expression.pdf"))
# boxplot of max 300 random samples
toplot <- correctedExpression[,sample(1:ncol(correctedExpression),min(ncol(correctedExpression), 300))]
boxplot(as.data.frame(toplot), main="expression of random samples", xlab="samples", ylab="expression")
hist(toplot, breaks=100)
# gene against gene correlation plots
cors <- cor(t(toplot))
hist(cors, breaks=100, main="correlation of genes", xlab="correlation")
dev.off()
