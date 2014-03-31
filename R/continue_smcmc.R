#########
### Method to continue MCMC estimation of a synMCMC object
#########

.continue.smcmc <- function(object, 
                            niter = object@niter,
                            nsim = object@nsim,
                            propCov = object@propCov, 
                            priorFun = object@priorFun,
                            targetRate = object@targetRate,
                            recompute = object@recompute,
                            multicore = object@multicore,
                            ncores = object@ncores,
                            cluster = NULL,
                            control = object@control,
                            ...)
{
  if(!is(object, "smcmc")) stop("To use mcmc you need an object of class \"smcmc\"")
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = object@control, outerCtrl = control)
  
  # For initPar and burn unless they have been specified by the user, we put
  # initPar to the final mcmc points in "object" and we don't do any more burn in.
  tmpObj <- smcmc(object = object, 
                  initPar = aaply(object@chains, 3, tail, n = 1),
                  niter = niter,
                  nsim = nsim,
                  propCov = propCov, 
                  burn = 0,
                  nchains = object@nchains,
                  priorFun = priorFun,
                  targetRate = targetRate,
                  recompute = recompute,
                  multicore = multicore,
                  ncores = ncores,
                  cluster = cluster,
                  control = ctrl,
                  ...)

  # Averaging old and new acceptance rate
  accRate <- (object@accRate * object@niter + tmpObj@accRate * tmpObj@niter) / (object@niter + tmpObj@niter)
  
  return(new(   "smcmc",
                tmpObj,
                
                initPar = object@initPar,   # Resetting the values of these two param, so we don't lose information  
                burn = as.integer(object@burn),
                
                niter = as.integer(object@niter + niter),
                
                accRate = accRate, 
                chains = abind(object@chains, tmpObj@chains, along = 1),
                llkChain = rbind(object@llkChain, tmpObj@llkChain),
                parStore = abind(object@parStore, tmpObj@parStore, along = 1),
                llkStore = rbind(object@llkStore, tmpObj@llkStore)
  ))
}

##############################################################
#' @details When \code{is("smcmc", object) == TRUE}  continues MCMC estimation of an object of class \code{smcmc}. All input parameters are defaulted to the corresponding
#' slots in the input object, with the exception of cluster. The chain restarts were it ended, burn-in is set to zero, the
#' same prior (if any) is used.
#' @aliases continue,smcmc-method
#' @rdname continue-generic
#' 
setMethod("continue", 
          signature = signature(object = "smcmc"), 
          definition = .continue.smcmc)



