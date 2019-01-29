# ------------------------------------------------------------------------------
#' Calculate peer factors for given data and covariates
#' Then get residual data matrix
#'
#' @param data nxg matrix (n=samples, g=genes/variables)
#' @param covariates nxc matrix (n=samples, c=covariates)
#' @param get.residuals Flag whether to directly return the residuals
#' calculated on the data matrix instead of the peer factors calculated
#' @param Nk Number of factors to estimate. Default: N/4
#'
#' @return Matrix of corrected expression data
#'
#' @author Johann Hawe
#'
#' @date 20170328
#'
# ------------------------------------------------------------------------------
correct_peer <- function(data, covariates=NULL,
                         Nk=ceiling(nrow(data)*0.25)) {

  if(!require(peer)) {
    stop("PEER needs to be installed to perform PEER normalization.")
  }

  # create model
  model <- PEER();

  # input has to be a matrix!
  PEER_setPhenoMean(model, as.matrix(data));

  # add the mean estimation as default since it is recommended in the tutorial
  # of peer. will return Nk+1 factors
  PEER_setAdd_mean(model, TRUE)

  # set number of hidden factors to identify. If unknown,
  # a good measure is N/4 (see howto)
  PEER_setNk(model, Nk);

  # should not be neccessary but increase anyways
  PEER_setNmax_iterations(model, 5000);

  if(!is.null(covariates)){
    # set matrix of known and important covariates,
    # since we want to acknowledge their effect
    PEER_setCovariates(model, as.matrix(covariates));
  }

  # learn
  PEER_update(model);

  re <- PEER_getResiduals(model)
  colnames(re) <- colnames(data)
  rownames(re) <- rownames(data)
  return(re)
}
