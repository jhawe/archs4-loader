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

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# -----------------------------------------------------------------
# Prepare libraries
# -----------------------------------------------------------------
# Check for dependencies and install if missing
print("Checking for required packages")
source("https://bioconductor.org/biocLite.R")

if(!require(rhdf5)) {
  biocLite("rhdf5")
  library(rhdf5)
}
if(!require(preprocessCore)) {
  biocLite("preprocessCore")
  library(preprocessCore)
}
if(!require(sva)) {
  biocLite("sva")
  library(sva)
}
if(!require(optparse)){
  biocLite("optparse")
  library(optparse)
}

# we assume these are here :P
library(ggplot2)
library(reshape2)

# -----------------------------------------------------------------
# Get snakemake params
# -----------------------------------------------------------------

# input
fsamples = snakemake@input$samples
fh5 = snakemake@input$h5

# output
fexpr = snakemake@output$expr
fplot = snakemake@output$plot
fdesign = snakemake@output$design

# params
# TODO: get these via snakemake
NORM = T
SVA = T
PEER = F
PLOT = F

# we normalize if we use PEER or SVA
NORM = NORM | PEER | SVA

if(SVA & PEER) {
  stop("Choose either COMBAT or PEER!")
}

if(PEER) {
  source("scripts/peer.R")
}

# -----------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------

# get the samples
source(fsamples)

# Retrieve information from compressed data
samples = h5read(fh5, "meta/Sample_geo_accession")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# get info for selected samples
tissue = h5read(fh5, "meta/Sample_source_name_ch1")
genes = h5read(fh5, "meta/genes")
series = h5read(fh5, "meta/Sample_series_id")
organism = h5read(fh5, "meta/Sample_organism_ch1")
molecule = h5read(fh5, "meta/Sample_molecule_ch1")
characteristics = h5read(fh5, "meta/Sample_characteristics_ch1")
description = h5read(fh5, "meta/Sample_description")
instrument = h5read(fh5, "meta/Sample_instrument_model")

design = cbind(sample=samples, tissue, series, organism, molecule, characteristics, description,
               instrument)
design = design[sample_locations,,drop=F]

# write design file
write.table(file=fdesign, design, sep="\t", quote=T, col.names=T, row.names=F)

# -----------------------------------------------------------------
# Extract and normalize expression data
# -----------------------------------------------------------------

# extract gene expression from compressed data
expression = h5read(fh5, "data/expression", index=list(1:length(genes), sample_locations))
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
  print("Removing batch effects using PEER.\nThis might take a while...")
  # transform the scaled counts to std normal per gene
  stdnorm <- function(x) {
    r = rank(x, ties.method="random")
    qnorm(r / (length(x) + 1))
  }
  transformed = apply(expression, 1, stdnorm)

  # peer expects an NxG matrix (N=#samples)
  # we use Nk=10 as in the lappaleinen et al. 2013 paper
  expression <- t(correct.peer(data=transformed, Nk=10))
  print("Done.")
}

# -----------------------------------------------------------------
# Save file and do some default plots
# -----------------------------------------------------------------

# Print file
write.table(expression, file=fexpr, sep="\t", quote=FALSE, col.names=NA, row.names=T)

if(PLOT) {

  # Create histogram of expression values and heatmap
  # of gene correlations
  print("Creating heatmap and expression histogram.")

  theme_set(theme_bw())

  pdf(fplot)
 
  # boxplot and histogram of max 150 random samples
  toplot <- data.frame(expression[,sample(1:ncol(expression),min(ncol(expression), 150))])
  toplot_melt <- melt(toplot)

  ggplot(toplot_melt, aes(y=value, x=variable)) + 
    geom_boxplot() + 
    xlab("samples") +
    ylab("normalized expression") + 
    ggtitle("Expression distribution of 150 random samples")

  ggplot(toplot_melt, aes(x=value)) + 
    geom_histogram(stat = "density") + 
    ggtitle("Distribution of expression values over 150 samples.")

  # gene against gene correlation plots
  corr <- cor(t(toplot))
  corr <- cbind.data.frame(correlation=corr[upper.tri(corr, diag=F)])
  ggplot(corr, aes(x=correlation)) + 
    geom_histogram() + 
    ggtitle("Distribution of gene-gene correlations.")

  dev.off()
}

# -----------------------------------------------------------------
# Session info
# -----------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
