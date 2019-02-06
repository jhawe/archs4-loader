# ------------------------------------------------------------------------------
#' Helper to quickly create a new design matrix from the h5 file
#'
#' @param fh5 The path to the h5 file from which to get the data
#' @param samp Subset of samples for which to create the file. Default: NULL
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
load_design <- function(fh5, samp=NULL) {

  # get info for all samples
  sample = h5read(fh5, "meta/Sample_geo_accession")
  tissue = h5read(fh5, "meta/Sample_source_name_ch1")
  genes = h5read(fh5, "meta/genes")
  series = h5read(fh5, "meta/Sample_series_id")
  organism = h5read(fh5, "meta/Sample_organism_ch1")
  molecule = h5read(fh5, "meta/Sample_molecule_ch1")
  characteristics = h5read(fh5, "meta/Sample_characteristics_ch1")
  description = h5read(fh5, "meta/Sample_description")
  instrument = h5read(fh5, "meta/Sample_instrument_model")
  processing = h5read(fh5, "meta/Sample_data_processing")
  library_selection = h5read(fh5, "meta/Sample_library_selection")
  library_source = h5read(fh5, "meta/Sample_library_source")
  library_strategy = h5read(fh5, "meta/Sample_library_strategy")
  platform_id = h5read(fh5, "meta/Sample_platform_id")
  relation = h5read(fh5, "meta/Sample_relation")
  status = h5read(fh5, "meta/Sample_status")
  title = h5read(fh5, "meta/Sample_title")
  taxid = h5read(fh5, "meta/Sample_taxid_ch1")
  type = h5read(fh5, "meta/Sample_type")
  
  # create design matrix
  design = cbind.data.frame(sample, tissue, series, organism, molecule,
                 characteristics, description,
                 instrument, processing, library_selection, library_source,
                 library_strategy, platform_id, relation, status, title,
                 taxid, type, typestringsAsFactors=F)
  
  # check whether we should extract only some specific samples
  if(!is.null(samp)) {
    design = dplyr::filter(design, sample %in% samp)
  }
  return(design)
}

# ------------------------------------------------------------------------------
#' Helper to get all samples in design which match the provided keywords
#'
#' @param design The path to the h5 file from which to get the data
#' @param keywords Vector of keywords to be matched. Can be 'NULL', in which
#' case all samples available in design will be returned. Default: NULL
#' @param exact Whether to exactly match the provided keywords or simply grep 
#' for them in the tissue anntotation of the design table. In the former case,
#' keyword results are combined with 'or', in the latter with 'and'.
#' Default: FALSE
#' @param filter_cancer Try to filter out cancer indicated samples? 
#' Default: FALSE
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
get_samples_from_design <- function(design, keywords=NULL, exact=FALSE,
                                    filter_cancer=FALSE) {

  # simply return all samples
  if(is.null(keywords)) {
    return(design$sample)
  }

  # define vector of which samples are to be used,
  # i.e. are matching the keywords
  use <- NULL
  for(k in keywords){
    if(is.null(use)){
      if(exact) {
        use <- grepl(paste0("^", k, "$"), design$tissue, ignore.case=F)
      } else {
        use <- grepl(k, design$tissue, ignore.case=T)
      }
    } else {
      if(exact) {
        use <- use |  grepl(paste0("^", k, "$"), design$tissue, ignore.case=F)
      } else {
        use <- use & grepl(k, design$tissue, ignore.case=T)
      }
    }
  }

  if(filter_cancer) {
    # try to filter out any samples which could be 'cancerous'
    use <- use & !grepl("cancer|carcino|adenom|blastom|tumour|tumor|sarcom", 
                         design$tissue, ignore.case=T)
    use <- use & !grepl("cancer|carcino|adenom|blastom|tumour|tumor|sarcom", 
                         design$description, ignore.case=T)
    use <- use & !grepl("cancer|carcino|adenom|blastom|tumour|tumor|sarcom", 
                         design$characteristics, ignore.case=T)
  }
  samples <- design[use,]$sample
  return(samples)
}

# ------------------------------------------------------------------------------
#' Process the raw expression data
#'
#' @param The expression matrix containing the raw counts
#' @param The normalization method to be applied. One of "quantile", "sva" or
#' "peer"
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
process_expression <- function(expression, norm_method) {

  # save colnames for resetting
  cnames <- colnames(expression)
  rnames <- rownames(expression)

  # normalize samples and correct for differences in gene count distribution
  expression = log2(expression+1)
  expression = normalize.quantiles(expression)

  # reset names after norm.quant
  rownames(expression) = rnames
  colnames(expression) = cnames

  # correct batch effects in gene expression
  # if we got some data...
  if(ncol(expression) > 10) {

    if(norm_method == "sva") {
      print("Removing batch effects using ComBat.")
      series = design$series[sample_locations]
      batchid = match(series, unique(series))
      expression <- ComBat(dat=expression, batch=batchid,
                           par.prior=TRUE, prior.plots=FALSE)
    } else if(norm_method == "peer") {
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
 return(expression)
}


# ------------------------------------------------------------------------------
#' Creates and saves t-SNE plot for raw and normalized expression matrices.
#'
#' @param raw The raw gene counts
#' @param norm The normalized gene expression matrix
#' @param tissues Tissue annotation for all samples in the matrices
#' @param fplot File to which to save the t-SNE plots
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
create_tsne <- function(raw, norm, tissues, fplot) {
  
  require(Rtsne)

  perp <- 30
  max_perp <- (nrow(raw) -1) / 3
  if(max_perp < perp){
    perp <- max_perp
    print(paste0("Using adjusted perplexity parameter: ", perp))
  }

  reduction <- Rtsne(raw, check_duplicates=FALSE, max_iter = 1000, theta = 0.0,
                     dims = 2, perplexity = perp)
  reduction2 <- Rtsne(norm, check_duplicates=FALSE, max_iter = 1000, theta = 0.0,
                     dims = 2, perplexity = perp)

  # get number of samples with respective gtex_tissues
  # and create new annotation for plotting
  counts <- table(tissues)
  tissues_wcounts <- paste(tissues, " (",
                           counts[tissues], ")",
                           sep="")


  pdf(fplot, width=15, height=12)

  toplot <- cbind.data.frame(reduction$Y, tissues_wcounts)
  colnames(toplot) <- c("dim1", "dim2", "label")
  plot_tsne(df=toplot,
            title="t-SNE on raw gene counts labeled \nby tissue information")

  toplot <- cbind.data.frame(reduction2$Y, tissues_wcounts)
  colnames(toplot) <- c("dim1", "dim2", "label")
  plot_tsne(df=toplot,
            title="t-SNE on normalized expression data labeled \nby tissue information")

  dev.off()

}
