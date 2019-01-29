# ------------------------------------------------------------------------------
#' The code has been adapted from the code given at
#' http://amp.pharm.mssm.edu/archs4/help.html
#'
#' We essentially extract the needed data from the downloaded
#' archive and normalize it if necessary. Additionally we create
#' some basic diagnostic plots on the data.
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
# Prepare libraries
# ------------------------------------------------------------------------------
library(rhdf5)
library(preprocessCore)
library(sva)
library(ggplot2)
library(reshape2)

# get helper methods
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# Get snakemake params
# ------------------------------------------------------------------------------

# input
fh5 = snakemake@input$h5

# output
fraw = snakemake@output$raw
fexpr = snakemake@output$expr
fplot = snakemake@output$plot
fdesign = snakemake@output$design

# params
keywords <- snakemake@wildcards$keywords

# check for special keyword
if("all_samples" %in% keywords) {
  # will skip the subsetting later
  keywords <- NULL
} else {
  # we expect "_" to be the word separator
  keywords <- strsplit(keywords, "_")[[1]]
}
# normalization params
norm_method <- snakemake@params$norm_method
if("sva" %in% norm_method) {
  SVA <- TRUE
  PEER <- FALSE
  NORM <- TRUE
} else if("peer" %in% norm_method) {
  PEER <- TRUE
  SVA <- FALSE
  NORM <- TRUE
} else if("quantile" %in% norm_method) {
  NORM <- TRUE
  PEER <- FALSE
  SVA <- FALSE
} else {
  # not recognized
  stop(paste0("Sorry, chosen normalization method '", norm_method,
              "' not supported."))
}

if(PEER) {
  source("scripts/peer.R")
}

# ------------------------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------------------------

# get samples
design <- load_design(fh5)
selected_samples <- get_samples_from_design(design, keywords)
print(paste0("Found ", length(selected_samples), " samples."))

# fail gracefully
if(length(selected_samples) < 1) {
  write.table(file=fdesign, design[c(),])
  write.table(file=fraw, 0)
  write.table(file=fexpr, 0)
  warning("No samples available for specified keywords.")
  quit("no", 0)
}

# write subsetted design file
design_subset <- dplyr::filter(design, sample %in% selected_samples)
write.table(file=fdesign, design_subset, sep="\t", quote=T, col.names=T, row.names=F)

# Get sample locations to subset expression data and the available genes
sample_locations = which(design$sample %in% selected_samples)

genes = h5read(fh5, "meta/genes")

# ------------------------------------------------------------------------------
# Extract and normalize expression data
# ------------------------------------------------------------------------------

# extract gene expression from compressed data
expression = h5read(fh5, "data/expression", 
                    index=list(1:length(genes), sample_locations))
H5close()

# set names
rownames(expression) = genes
colnames(expression) = design$sample[sample_locations]

# write raw expression
write.table(expression, file=fraw, sep="\t",
            quote=FALSE, col.names=NA, row.names=T)

# normalize samples and correct for differences in gene count distribution
if(NORM | SVA | PEER){
  print("Normalizing data.")
  expression = log2(expression+1)
  expression = normalize.quantiles(expression)
}

# reset names after norm.quant
rownames(expression) = genes
colnames(expression) = design$sample[sample_locations]

# correct batch effects in gene expression
# if we got some data...
if(ncol(expression) > 10) {

  if(SVA) {
    print("Removing batch effects using ComBat.")
    series = design$series[sample_locations]
    batchid = match(series, unique(series))
    expression <- ComBat(dat=expression, batch=batchid, 
                         par.prior=TRUE, prior.plots=FALSE)
  } else if(PEER) {
    print("Removing batch effects using PEER.")
    # transform the scaled counts to std normal per gene
    stdnorm <- function(x) {
      r = rank(x, ties.method="random")
      qnorm(r / (length(x) + 1))
    }
    transformed = apply(expression, 1, stdnorm)

    # peer expects an NxG matrix (N=#samples)
    # we use Nk=10 as in the lappaleinen et al. 2013 paper
    expression <- t(correct_peer(data=transformed, Nk=10))
    print("Done.")
  }

}
# ------------------------------------------------------------------------------
# Save file
# ------------------------------------------------------------------------------
write.table(expression, file=fexpr, sep="\t",
            quote=FALSE, col.names=NA, row.names=T)

# ------------------------------------------------------------------------------
# Session info
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
