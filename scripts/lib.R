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
  samples = h5read(fh5, "meta/Sample_geo_accession")
  tissue = h5read(fh5, "meta/Sample_source_name_ch1")
  genes = h5read(fh5, "meta/genes")
  series = h5read(fh5, "meta/Sample_series_id")
  organism = h5read(fh5, "meta/Sample_organism_ch1")
  molecule = h5read(fh5, "meta/Sample_molecule_ch1")
  characteristics = h5read(fh5, "meta/Sample_characteristics_ch1")
  description = h5read(fh5, "meta/Sample_description")
  instrument = h5read(fh5, "meta/Sample_instrument_model")

  # create design matrix
  design = cbind.data.frame(sample=samples, tissue, series, organism, molecule,
                 characteristics, description,
                 instrument, stringsAsFactors=F)
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
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
get_samples_from_design <- function(design, keywords=NULL, exact=FALSE) {

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
        use <- grepl(paste0("^", k, "$"), design$tissue, ignore.case=T)
      } else {
        use <- grepl(k, design$tissue, ignore.case=T)
      }
    } else {
      if(exact) {
        use <- use |  grepl(paste0("^", k, "$"), design$tissue, ignore.case=T)
      } else {
        use <- use & grepl(k, design$tissue, ignore.case=T)
      }
    }
  }

  samples <- design[use,]$sample
  return(samples)
}

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
