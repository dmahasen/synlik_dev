########
#' Extract estimates from \code{mcmc} objects
#'
#' @param object  ("mcmc") object.
#' @param method  if set to "emp" then the output will be the empirical covariance matrix
#'                of the chain(s). In set to "asym" then a quadratic model will be fitted
#'                through the likelihood evaluations and the output covariance will be the 
#'                inverse of the negative Hessian.  
#' @param nreps   (integer) number of simulations used to tilt the negative Hessian (-H)
#'                toward positive definiteness (PD). Used only if the -H is not PD and \code{method == "asym"}.
#' @param boot    (logical) relevant only if -H is not PD. 
#'                If TRUE hessians will be simulated by resampling parameters and likelihoods. 
#'                If FALSE hessians will be simulated the asymptotic distribution of the regression
#'                coefficients.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @method coef mcmc
#' @export
#' 

coef.mcmc <- function(object,
                      method = "bayes",
                      
                      burn = 0,
                      quant = 0.1,
                      
                      nreps = 1000, 
                      boot = TRUE, 
                      ...)
{
  if( !(method %in% c("bayes", "mle")) ) stop("method should be either \"mle\" or \"bayes\"")
  
  # Discarding initial part of the chains
  if( burn > 0 ){
    object$llkStore <- object$llkStore[-(1:burn), , drop = FALSE]
    object$parStore <- object$parStore[-(1:burn), ,  , drop = FALSE]
  }
  
  if(method == "mle")
  {
    
    varPar <- diag(object$propCov) > 0
    parEstim <- object$initPar[1, ]
    
    regr <- .quadLikWeighRegr(llk = object$llkStore, 
                              parMat = object$parStore[ , varPar, , drop = FALSE],
                              nreps = nreps, 
                              boot = boot,
                              quant = quant,
                              ...)
    
    # Maximizing [ f = a + t(b) * x + t(x) * C * x ] with respect to x
    parEstim[varPar] <- - 0.5 * solve(regr$quad, regr$linear[-1])
     
  } else {
    
    parEstim <- colMeans( aaply(object$chains, 3, colMeans, .drop = FALSE) )
    
  }
  
  return(parEstim) 
}

