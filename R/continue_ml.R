#########
### Method to continue ml() estimation of an ml object
#########

# For initPar, initCov amd control unless they have been specified by the user, we use the final 
# values contained in object (ex: initPar is the final value in "object")
continue.ml <- function(object, 
                        initCov = object$initCov,
                        np = object$np,
                        niter = object$niter,
                        alpha = object$alpha,
                        temper = rep(1, niter),
                        recycle = object$recycle,
                        multicore = object$multicore,
                        ncores = object$ncores,
                        cluster = NULL,
                        constr = object$constr,
                        verbose = object$verbose,
                        ...)
{
    
  stop("this function cannot work at the moment")
  
  outObj <- ml(likFun = object$likFun,
               initPar = drop( tail(object$estim, 1) ),
               initCov = initCov * (object$alpha ^ object$niter),
               np = np,
               niter = niter,
               priorFun = object$priorFun,
               alpha = alpha,
               temper = temper,
               recycle = recycle,
               multicore = multicore,
               ncores = ncores,
               constr = object$constr,
               cluster = cluster,
               ...)
  
  outObj$initPar <- object$initPar
  outObj$initCov <- object$initCov
  outObj$niter <- outObj$niter + object$initCov
  
  outObj$estim <- rbind(object$estim, outObj$estim)
  outObj$simLogLik <- append(object$simLogLik, outObj$simLogLik)
  outObj$simLogPrior <- append(object$simLogPrior, outObj$simLogPrior)
  outObj$simPar <- rbind(object$simPar, outObj$simPar)
  
  return( outObj )
}