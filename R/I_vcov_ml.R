########
#' Parameters covariance matrix for "sml" objects
#'
#' @param llk     (numeric) vector containing the log-likelihood values.
#' @param parMat  (matrix) where the i-th row is the parameter vector corresponding to the i-th log-likelihood.
#' @param nreps   (integer) number of simulations used to tilt the negative Hessian (-H)
#'                toward positive definiteness. Used only if -H is not PD.
#' @param boot    (logical) relevant only if -H is not PD. 
#'                If TRUE hessians will be simulated by resampling parameters and likelihoods. 
#'                If FALSE hessians will be simulated the asymptotic distribution of the regression
#'                coefficients.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' 

.vcov.ml <- function(llk, parMat, nreps = 1000, boot = TRUE, burn = 0, quant = 0.1, ...)
{ 

  # Discarding initial part of the chains
  if( burn > 0 ){
    llk <- llk[-(1:burn)]
    parMat <- parMat[-(1:burn), , drop = FALSE]
  }
  
  hess <- .quadLikWeighRegr(llk = llk, 
                            parMat = parMat, 
                            nreps = nreps, 
                            boot = boot,
                            quant = quant,
                            ...)$hessian 
  
  # Getting covariance, standard errors and confidence intervals
  covar <- .qrInverse(-hess)
  rownames(covar) <- colnames(covar) <- colnames(parMat)
  
  return(covar) 
}

