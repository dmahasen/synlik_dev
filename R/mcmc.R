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
  
  # Preparing input for the chains: each core get a list where the first element is a vector
  # with the initial values, and the second is the initial value of the log-likelihood
  if( is.null(control$initLoglik) ) control$initLoglik <- rep(-1e99, nchains)
  initSettings <- mapply("list", split(initPar, 1:nchains), control$initLoglik, SIMPLIFY = FALSE)

  # multicore, ncores and cluster go in the "..."
  singleChain <- function(input, ...)
  {
    control$initLoglik <- input[[2]]
    
    .mcmc(likFun = likFun, 
          initPar = input[[1]], 
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
    out <- llply(.data = initSettings,
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
  
  return( list("likFun"  = likFun,
               "initPar" = initPar,
               "niter"   = niter,
               "propCov" = propCov,
               "burn"    = burn,
               "nchains" = nchains,
               "priorFun" = priorFun,
               "targetRate" = targetRate,
               "recompute" = recompute,
               "multicore" = multicore,
               "ncores" = ncores,
               "control" = control,
               
               "accRate" = accRate, 
               "chains" = chains, 
               "llkChain" = llkChain,
               "parStore" = parStore,
               "llkStore" = llkStore) )
}


##############
## INTERNAL
##############

.mcmc <- function(likFun, 
                  initPar, 
                  niter, 
                  propCov, 
                  burn = 0,
                  priorFun = function(param, ...) 0,
                  targetRate = NULL,
                  recompute = FALSE,
                  multicore = !is.null(cluster),
                  cluster = NULL,
                  ncores = detectCores() - 1, 
                  control = list(),
                  ...)
{
  
  totalIter <- niter + burn
  
  # Control list which will be used internally
  ctrl <- list( "theta" = 0.5,
                "adaptStart" = 0,
                "adaptStop" = totalIter,
                "saveFile" = NULL,
                "saveFreq" = 100,
                "verbFreq" = 500, 
                "verbose"  = FALSE,
                "initLoglik" = -1e99)
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  # Safety checks on ctrl
  if(ctrl$theta > 1 || ctrl$theta < 0.5) stop("control$theta should be between 0.5 and 1")
  if( !is.null(ctrl$saveFile) && !is.character(ctrl$saveFile) ) stop("\"ctrl$saveFile\" should be a character vector")
  stopifnot( ctrl$adaptStart <= ctrl$adaptStop, 
             ctrl$adaptStart >= 0,  
             ctrl$adaptStop <= totalIter) 
  
  # Check other arguments
  if(is.matrix(propCov) == FALSE) propCov <- diag(propCov)
  if(nrow(propCov) != ncol(propCov)) stop("propCov should be a square matrix")
  if(nrow(propCov) != length(initPar)) stop("nrow(propCov) != length(initPar)")
  
  # If a parameter has variance 0 in the proposal we save it's index in "fixPar"
  # we modidy the covariance and we save the initial covariance
  fixPar <- which( diag(propCov) == 0 )
  anyFix <- ( length(fixPar) > 0 )
  savedCov <- propCov
  if(anyFix) diag(propCov)[fixPar] <- 0.1
  
  cholFact <- t( chol( unname(propCov) ) )
  
  nPar <- length(initPar);
  currPar <- unname( initPar );
  propPar <- numeric(nPar);
  
  mcmcSample <- parStore <- matrix(NA, niter, nPar);
  llkChain <- llkStore <- numeric(niter);
  
  currLogLik <- propLogLik <- tmpLik <- ctrl$initLoglik;
  currPrior <- priorFun(initPar, ...)
  
  accept <- 0
  
  unifVar <- runif(totalIter)
  
  # Mcmc main loop
  istore <- 1
  for(ii in 1:totalIter){
    
    # Propose new parameters
    pert <- rnorm(nPar)
    propPar <- currPar + as.vector( cholFact %*% pert )
    
    # Fix some parameters (if necessary) and check prior
    if( anyFix ) propPar[fixPar] <- currPar[fixPar]
    propPrior <- priorFun(propPar, ...)
    
    if( is.finite(propPrior) )
    { 
      # Compute likelihood of proposed param
      propLogLik <- try( likFun(param = propPar, multicore = multicore, ncores = ncores, cluster = cluster, ...) )
      if( !is.numeric(propLogLik) || !is.finite(propLogLik) ) {
        warning(paste("One function was equal to", propLogLik, "and I put it to -Inf."))
        propLogLik <- -Inf 
      }
      
      # (Optionally) recompute likelihood at old parameters
      if(recompute){ 
        tmpLik <- try( likFun(param = currPar, multicore = multicore, ncores = ncores, cluster = cluster, ...) )
        if( is.numeric(tmpLik) && is.finite(tmpLik) ){ 
          currLogLik <- tmpLik 
        } else {
          warning(paste("One function was equal to", currLogLik, "and I put it to -Inf."))
        }
      }
      
      # Compute acceptance ratio
      likDiff <- propLogLik  - currLogLik + propPrior - currPrior
      alpha <- min(1, exp(likDiff))
      if( !is.finite(alpha) ) alpha <- ifelse( likDiff >= 0, 1, 0) 
      
      # Accept/Reject
      if ( likDiff > log(unifVar[ii]) ) {
        currPar <- propPar
        currPrior <- propPrior
        currLogLik <- propLogLik
        if(ii > burn) accept <- accept + 1
      }
      
    } else { alpha <- 0 }
    
    # Store iteration if iteration > burn-in
    if(ii > burn) {
      mcmcSample[istore, ] <- currPar;
      llkChain[istore] <- currLogLik;
      parStore[istore, ] <- propPar
      llkStore[istore] <- propLogLik
      istore <- istore + 1
    }
    
    # (Optionally) adapt the proposal distribution, by updatint the transpose of its Cholesky factor
    if( !is.null(targetRate) && (ii >= ctrl$adaptStart) && (ii <= ctrl$adaptStop) )
    {
      cholFact <- .adaptChol(nPar = nPar, iter = ii, S = cholFact, U = pert, 
                             gamma = ctrl$theta, alpha = alpha, accRate = targetRate)
    }
    
    # (Optionally) save the object to file
    if( !is.null(ctrl$saveFile) && !(ii %% ctrl$saveFreq) ){ 
      
      out <- list("accRate" = accept / ii, "chains"  = mcmcSample, "llkChain" = llkChain)
      
      save(file = ctrl$saveFile, out)   
    }
    
    # (Optionally) print out intermediate results
    if( ctrl$verbose == TRUE && (ii > burn) && !(ii %% ctrl$verbFreq) )
    {
      tmp <- colMeans(mcmcSample[1:ii, ], na.rm = TRUE)
      names(tmp) <- names(object@param)
      cat(paste("Empirical posterior means at iteration", ii - burn, "\n"))
      print(tmp)
    }
    
  }
  
  colnames(mcmcSample) <- names(initPar)
  
  return( list("accRate" = accept / niter, 
               "chains"  = mcmcSample, 
               "llkChain" = llkChain, 
               "parStore" = parStore,
               "llkStore" = llkStore) )
  
}
