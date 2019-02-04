# ------------------------------------------------------------------------------
#' The code has been adapted from the code given at
#' http://amp.pharm.mssm.edu/archs4/help.html
#'
#' We essentially extract the needed data from the downloaded
#' archive and normalize it if necessary. Additionally we create
#' some basic diagnostic plots on the data.
#'
#' This script differs from the 'extract_data.R' script in that it expacts
#' keywords to be passed via the snakemake parameters and all these keywords
#' are matched exactly (rather than simply grepped) in the tissue meta data 
#' annotation.
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
library(Rtsne)

# get helper methods
source("scripts/lib.R")
source("scripts/plotting.R")

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
keywords <- snakemake@params$keywords
# we expect "_" to be the word separator
keywords <- strsplit(keywords, "\\|")[[1]]

# ------------------------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------------------------

# get samples
design <- load_design(fh5)
selected_samples <- get_samples_from_design(design, keywords, exact=T)
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

expression_processed <- process_expression(expression, "sva")

# ------------------------------------------------------------------------------
print("Saving expression data.")
# ------------------------------------------------------------------------------
write.table(expression_processed, file=fexpr, sep="\t",
            quote=FALSE, col.names=NA, row.names=T)

# ------------------------------------------------------------------------------
print("Creating tSNE and plotting.")
# ------------------------------------------------------------------------------

# ensure column sort
raw <- t(expression[,design_subset$sample])
norm<- t(expression_processed[,design_subset$sample])

reduction <- Rtsne(raw, check_duplicates=FALSE, max_iter = 1000, theta = 0.0,
                   dims = 2, perplexity = 30)
reduction2 <- Rtsne(norm, check_duplicates=FALSE, max_iter = 1000, theta = 0.0,
                   dims = 2, perplexity = 30)

# get number of samples with respective gtex_tissues
# and create new annotation for plotting
counts <- table(design_subset$tissue)
tissues_wcounts <- paste(design_subset$tissue, " (",
                         counts[design_subset$tissue], ")",
                         sep="")

  ggtitle("t-SNE on raw gene expression data labeled \nby gtex_tissue information")

pdf(fplot, width=15, height=12)
y <- reduction$Y
plot_tsne(y[,1], y[,2], tissues_wcounts, 
          "t-SNE on raw gene counts labeled \nby tissue information")
y <- reduction2$Y
plot_tsne(y[,1], y[,2], tissues_wcounts, 
          "t-SNE on normalized expression data labeled \nby tissue information")
dev.off()

# ------------------------------------------------------------------------------
# Session info
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
