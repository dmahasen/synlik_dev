#######
### Method to estimate the parameters through MCMC
#######
# General mcmc procedure.
# Arguments that differ from smcmc:
# object (missing)
# nsim   (missing)
# libraries = names of libraries that will be loaded into the cluster
mcmc <- function(likFun, 
                 initPar, 
                 niter, 
                 propCov, 
                 burn = 0,
                 nchains = 1,
                 priorFun = function(param, ...) 0,
                 targetRate = NULL,
                 recompute = FALSE,
                 multicore = !is.null(cluster),
                 cluster = NULL,
                 ncores = detectCores() - 1, 
                 libraries = c(),
                 control = list(),
                 ...)
{
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available to funToApply on cluster
    .forceEval(ALL = TRUE)
    
      tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = libraries, exportALL = TRUE)
      cluster <- tmp$cluster
      ncores <- tmp$ncores
      clusterCreated <- tmp$clusterCreated
  }
  
  # multicore, ncores and cluster go in the ...
  funToApply <- function(notUsed, ...)
  {
    .mcmc(likFun = likFun, 
          initPar = initPar, 
          niter = niter, 
          propCov = propCov, 
          burn = burn,
          priorFun = priorFun,
          targetRate = targetRate,
          recompute = recompute,
          control = control,
          ...) 
  }
  
  out <-  if(nchains > 1 && multicore)
  {
    # nchains each on one cores, each node of the cluster take one chain
    parLapply(cl  = cluster,
              X   = 1:nchains, 
              fun = funToApply,
              multicore = FALSE,
              ncores = 1,
              cluster = NULL,
              ...)
  } else {
    # Only one chain which uses the whole cluster
    lapply(1:nchains, 
           funToApply,
           multicore = multicore,
           ncores = ncores,
           cluster = cluster,
           ...)
  }
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  if(nchains == 1) out <- out[[1]]
  
  return(out)
}
