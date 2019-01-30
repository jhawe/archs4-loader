# ------------------------------------------------------------------------------
#' The code has been adapted from the code given at
#' http://amp.pharm.mssm.edu/archs4/help.html
#'
#' Here we get a set of gtex tissues for which we can directly find samples
#' in archs4 and investigate the tSNE again
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
ftissues = snakemake@tissues
tissues <- read.table(ftissues, header=F, stringsAsFactors=F)[,1]

# output
fraw = snakemake@output$raw
fexpr = snakemake@output$expr
fplot = snakemake@output$plot
fdesign = snakemake@output$design

# params
norm_method <- snakemake@params$norm_method
if(!norm_method %in% c("sva", "peer", "quantile")) {
  # not recognized
  stop(paste0("Sorry, chosen normalization method '", norm_method,
              "' not supported."))
}

if(norm_method == "peer") {
  source("scripts/peer.R")
}

# ------------------------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------------------------

# get samples
design <- load_design(fh5)
selected_samples <- c()
for(tissue in tissues) {
  tis <- strsplit(tissue, "_")[[1]]
  samp <- get_samples_from_design(design, tis)
  names(samp) <- tissue
  selected_samples <- c(selected_samples, samp)
}
print(paste0("Found ", length(selected_samples), " samples."))

# write subsetted design file
design_subset <- dplyr::filter(design, sample %in% selected_samples) %>%
  dplyr::mutate(gtex_tissue = names(selected_samples))

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

expression_processed <- process_expression(expression, norm_method)

# directly perform tSNE and save the plot
# ensure column sort
raw <- dplyr::select(raw, one_of(design_subset$sample))

reduction <- Rtsne(raw, check_duplicates=FALSE, max_iter = 1000, theta = 0.0,
                   dims = 2, perplexity = perp)

# plotting
toplot <- reduction$Y
colnames(toplot) <- c("dim1", "dim2")

toplot <- cbind(toplot, tissue = design_subset$tissue,
                instrument = design_subset$instrument,
                series=design_subset$series,
                gtex_tissue = design_subset$gtex_tissue)

gp1 <-ggplot(toplot, aes(x=dim1, y=dim2, col=tissue)) +
  geom_point() +
  facet_wrap(~type, nrow=2) +
  ggtitle("t-SNE on gene expression data labeled \nby 'tissue' meta-data information")

pdf(fplot)
gp1
dev.off()

# ------------------------------------------------------------------------------
# Save file
# ------------------------------------------------------------------------------
write.table(expression_processed, file=fexpr, sep="\t",
            quote=FALSE, col.names=NA, row.names=T)

# ------------------------------------------------------------------------------
# Session info
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
