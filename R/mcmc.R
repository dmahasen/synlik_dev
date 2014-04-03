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
                 control = list(),
                 ...)
{
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available to singleChain on cluster
    .forceEval(ALL = TRUE)
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  # Prepare initial values in matrix form for apply(). Each row is a parameters set.
  if( is.vector(initPar) ) initPar <- t(initPar)
    
  if( is.matrix(initPar) ){
    
    if( nrow(initPar) == 1 ) initPar <- matrix(initPar, nchains, ncol(initPar), 
                                               byrow = TRUE, dimnames = list(NULL, colnames(initPar)))
    
    if( nrow(initPar) != nchains ) stop("nrow(initPar) should be either 1 or nchains")
    
  } else {
    
    stop("initPar should be either a matrix or a vector")
         
  }
    
  # multicore, ncores and cluster go in the ...
  singleChain <- function(input, ...)
  {
    .mcmc(likFun = likFun, 
          initPar = input, 
          niter = niter, 
          propCov = propCov, 
          burn = burn,
          priorFun = priorFun,
          targetRate = targetRate,
          recompute = recompute,
          control = control,
          ...) 
  }
  
  # Each chain goes on one node
  if(nchains > 1 && multicore){
    funMulti <- FALSE
    funCluster <- FALSE
    funNcores <- 1
  } else {
  # Only one chains where likelihood is evaluated in parallel 
    funMulti <- multicore
    funCluster <- cluster
    funNcores <- ncores
  }
  
  # Launch MCMC chain(s)
  withCallingHandlers({
    out <- alply(.data = initPar,
                 .margins = 1,
                 .fun = singleChain,
                 .parallel = multicore && (nchains > 1),
                 # ... from here
                 multicore = funMulti,
                 ncores = funNcores,
                 cluster = funCluster,
                 ...
                 )
  }, warning = function(w) {
    # There is a bug in plyr concerning a useless warning about "..."
    if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
      invokeRestart("muffleWarning")
  })
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  # Extracting acceptance rates
  accRate <- sapply(out, "[[", "accRate")
  
  # Putting chains in an niter X npar X nchains array
  chains <- lapply(out, "[[", "chains")
  chains <- do.call("abind", c(chains, "along" = 3))
  
  # Putting chains in an niter X npar X nchains array
  parStore <- lapply(out, "[[", "parStore")
  parStore <- do.call("abind", c(parStore, "along" = 3))
  
  # Putting llkChain is a matrix niter X nchains
  llkChain <- do.call("cbind", lapply(out, "[[", "llkChain")) 
  llkStore <- do.call("cbind", lapply(out, "[[", "llkStore"))
  
  return( list("initPar" = initPar,
               "propCov" = propCov,
                
               "accRate" = accRate, 
               "chains" = chains, 
               "llkChain" = llkChain,
               "parStore" = parStore,
               "llkStore" = llkStore) )
}
