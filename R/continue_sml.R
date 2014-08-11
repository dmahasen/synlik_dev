

.continue.sml <- function(object, 
                          initCov = object@initCov, 
                          np = object@np, 
                          nsim = object@nsim, 
                          niter = object@niter, 
                          alpha = object@alpha, 
                          priorFun = object@priorFun,
                          temper = rep(1, niter), 
                          recycle = object@recycle,
                          multicore = object@multicore, 
                          ncores = object@ncores, 
                          cluster = NULL, 
                          constr = object@constr, 
                          verbose = FALSE, 
                          ...)
{
  # Force evaluation of everything in the environment, so it will available to likfun on cluster
  if( multicore ) .forceEval(ALL = TRUE)
  
  # Function that will be used by sapply() or clusterApply to evaluate the likelihood
  likFun <- function(param, ...)
  {
    require("synlik")
    slik(object, param, nsim, multicore = FALSE, cluster = NULL, ...)
  }
  
  if( identical(initCov, object@initCov) ) initCov <- initCov * (object@alpha ^ object@niter)
  
  # Calling general maximum likelihood method
  tmp <- ml(likFun = likFun,
            initPar = drop( tail(object@estim, 1) ),
            initCov = initCov,
            np = np,
            niter = niter,
            alpha = alpha,
            priorFun = priorFun,
            temper = temper,
            recycle = recycle,
            multicore = multicore,
            ncores = ncores,
            constr = constr,
            cluster = cluster,
            ...)
  
  colnames(tmp$estim) <- names(object@initPar)
  
  return( .sml(object, 
               initPar = object@initPar, 
               initCov = object@initCov, 
               np = as.integer(np), 
               nsim = as.integer(nsim), 
               niter = as.integer(object@niter + niter),
               priorFun = priorFun,
               alpha = alpha,
               temper = temper,
               recycle = recycle,
               multicore = multicore, 
               ncores = as.integer(ncores), 
               constr = constr,
               
               estim = rbind(object@estim, tmp$estim),
               simLogLik = append(object@simLogLik, tmp$simLogLik),
               simLogPrior = append(object@simLogPrior, tmp$simLogPrior),
               simPar = rbind(object@simPar, tmp$simPar))  )
}


##############################################################
#' @details When \code{is("sml", object) == TRUE}  continues SML estimation of an object of class \code{sml}. All input parameters are defaulted to the corresponding
#' slots in the input object, with the exception of cluster. The optimizer restarts were it ended.
#' @aliases continue,sml-method
#' @rdname continue-generic
#' 
setMethod("continue", 
          signature = signature(object = "sml"), 
          definition = .continue.sml)
