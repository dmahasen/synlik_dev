
sml <- function(object, initPar, initCov, np, nsim, niter, alpha = 0.95,
                priorFun = NULL, temper = NULL, recycle = FALSE,
                multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                constr = list(), verbose = FALSE, ...)
{
  if( is.null(names(initPar)) ) names(initPar) <- names(object@param)
  
  # Force evaluation of everything in the environment, so it will available to likfun on cluster
  if( multicore ) .forceEval(ALL = TRUE)
  
  # Function that will be used by sapply() or clusterApply to evaluate the likelihood
  likFun <- function(param, ...)
  {
    require("synlik")
    slik(object, param, nsim, multicore = FALSE, cluster = NULL, ...)
  }
    
  # Calling general maximum likelihood method
  tmp <- ml(likFun = likFun, 
            initPar = initPar, 
            initCov = initCov, 
            np = np, 
            niter = niter,
            priorFun = priorFun,
            alpha = alpha,
            temper = temper,
            recycle = recycle,
            multicore = multicore, 
            ncores = ncores, 
            cluster = cluster, 
            constr = constr,
            verbose = verbose,
            ...)
  
  colnames(tmp$estim) <- names(initPar)
  
  return( .sml(object, 
               initPar = initPar, 
               initCov = initCov, 
               np = as.integer(np), 
               nsim = as.integer(nsim), 
               niter = as.integer(niter),
               priorFun = priorFun,
               alpha = alpha,
               temper = temper,
               recycle = recycle,
               multicore = multicore, 
               ncores = as.integer(ncores), 
               constr = constr,
               
               estim = tmp$estim,
               simLogLik = tmp$simLogLik,
               simLogPrior = tmp$simLogPrior,
               simPar = tmp$simPar)  )
}









########
# Adaptation
########
#     if( adapt )
#     {
#       propUpdate <- matrix(0, length(parMean), length(parMean))
#       for(kk in 1:np)
#       {
#         propUpdate <- propUpdate + w[kk] * tcrossprod(simPar[kk, ] - parMean, simPar[kk, ] - parMean)
#       }
#       
#       ESS <- 1 / sum(w ^ 2)
#       beta <- min(ESS / length(parMean), 0.25)
#       
#       parCov <- alpha * ( (1 - beta) * parCov + beta * propUpdate )
#     }