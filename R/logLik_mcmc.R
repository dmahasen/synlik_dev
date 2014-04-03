########
#' Extract maximized log-likelihood for \code{mcmc} objects
#'
#' @param object  ("mcmc") object.
#' @param nreps   (integer) number of simulations used to tilt the negative Hessian (-H)
#'                toward positive definiteness (PD). Used only if the -H is not PD and \code{method == "asym"}.
#' @param boot    (logical) relevant only if -H is not PD. 
#'                If TRUE hessians will be simulated by resampling parameters and likelihoods. 
#'                If FALSE hessians will be simulated the asymptotic distribution of the regression
#'                coefficients.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @method logLik mcmc
#' @export
#' 

logLik.mcmc <- function(object,
                        
                        burn = 0,
                        quant = 0.1,
                        
                        nreps = 1000, 
                        boot = TRUE, 
                        ...)
{
  
  # Discarding initial part of the chains
  if( burn > 0 ){
    object$llkStore <- object$llkStore[-(1:burn), , drop = FALSE]
    object$parStore <- object$parStore[-(1:burn), varPar,  , drop = FALSE]
  }
  
  varPar <- diag(object$propCov) > 0
  parMax <- object$initPar[1, ]
  
  regr <- .quadLikWeighRegr(llk = object$llkStore, 
                            parMat = object$parStore[ , varPar, , drop = FALSE],
                            nreps = nreps, 
                            boot = boot,
                            quant = quant,
                            ...)
  
  # Maximizing [ f = a + t(b) * x + t(x) * C * x ] with respect to x
  parMax[varPar] <- - 0.5 * solve(regr$quad, regr$linear[-1])
    
  llkMax <- regr$linear[1] + t(regr$linear[-1]) %*% parMax[varPar] + t(parMax[varPar]) %*% regr$quad %*% parMax[varPar]
  
  attr(llkMax, "df") <- length(varPar)
  class(llkMax) <- "logLik"
  
  return(llkMax) 
}
