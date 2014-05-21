#########
### Method to continue MCMC estimation of a "mcmc" object
#########

continue.mcmc <- function(object, 
                          niter = object$niter,
                          propCov = object$propCov, 
                          priorFun = object$priorFun,
                          targetRate = object$targetRate,
                          recompute = object$recompute,
                          multicore = object$multicore,
                          ncores = object$ncores,
                          cluster = NULL,
                          control = object$control,
                          ...)
{
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = object$control, outerCtrl = control)
  
  # The initial log-likelihood has to be equal to the last one in the previous part of the chain
  ctrl$initLoglik <- drop( tail(object$llkChain, 1) )
  
  # For initPar and burn unless they have been specified by the user, we put
  # initPar to the final mcmc points in "object" and we don't do any more burn in.
  tmpObj <- mcmc(likFun = object$likFun, 
                 initPar = aaply(object$chains, 3, tail, n = 1),
                 niter = niter,
                 propCov = propCov, 
                 burn = 0,
                 nchains = object$nchains,
                 priorFun = priorFun,
                 targetRate = targetRate,
                 recompute = recompute,
                 multicore = multicore,
                 ncores = ncores,
                 cluster = cluster,
                 control = ctrl,
                 ...)
  
  # Averaging old and new acceptance rate
  accRate <- (object$accRate * object$niter + tmpObj$accRate * tmpObj$niter) / (object$niter + tmpObj$niter)
  
  return(list(  "likFun" = object$likFun,
                
                "initPar" = object$initPar,   # Resetting the values of these two param, so we don't lose information  
                "burn" = as.integer(object$burn),
                "niter" = as.integer(object$niter + niter),
                "nchains" = object$nchains,
                "propCov" = propCov,
                "priorFun" = priorFun,
                "targetRate" = targetRate,
                "recompute" = recompute,
                "multicore" = multicore,
                "ncores" = ncores,
                
                "accRate" = accRate, 
                "chains" = abind(object$chains, tmpObj$chains, along = 1),
                "llkChain" = rbind(object$llkChain, tmpObj$llkChain),
                "parStore" = abind(object$parStore, tmpObj$parStore, along = 1),
                "llkStore" = rbind(object$llkStore, tmpObj$llkStore)
  ))
}
