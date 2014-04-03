

vcov.mcmc <- function(object,
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
    object$parStore <- object$parStore[-(1:burn), varPar,  , drop = FALSE]
  }
  
  if(method == "mle")
  {
    varPar <- diag(object$propCov) > 0
    covar <- object$propCov * 0
    hess <-  .quadLikWeighRegr(llk = object$llkStore, 
                               parMat = object$parStore[ , varPar, , drop = FALSE], 
                               nreps = nreps, 
                               boot = boot, 
                               quant = quant,
                               ...)$hessian 
    
    covar[varPar, varPar] <- .qrInverse(-hess)
    
  } else {
    
    covar <- alply(object$parStore, 3, cov)
    covar <- Reduce("+", covar) / length(covar)
    
  }
  
  covar[ is.na(covar) ] <- 0
  
  rownames(covar) <- colnames(covar) <- names(object$initPar)
  
  return(covar) 
}