
# Selecting initial covariance for models of class SynMaxlik, visually

selectCov <- function(object, 
                      ranges,
                      param = object@param,
                      nreps, nsim, 
                      nBoot = 0,
                      addRegr = TRUE, 
                      constr = list(), 
                      multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                      verbose = FALSE,
                      ...)
{
  
  # Save the simulator, the constraints and the parameters
  savedSimulator <- object@simulator 
  force(savedSimulator)
  
  fixConstr <- constr
  fixParam <- param
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available to singleChain on cluster
    .forceEval(ALL = TRUE)
    verbose <- FALSE
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  grads <- list()
  
  # Loop trough the parameters for which a range of standard deviation has been defined
  for(parNam in names(ranges))
  {
    # Index of the parameter currently considered
    parIndex <- which( names(object@param) == parNam )
    
    # Create a simulator that takes only one parameter has input
    object@simulator <- function(param, nsim, extraArgs, ...){ 
      
      parMat <- matrix(fixParam, nrow(param), length(fixParam), byrow = TRUE)
      
      parMat[ , parIndex] <- param 
      
      savedSimulator(parMat, nsim, extraArgs, ...)
      
    }
    
    # Using only the constraint of the current parameter
    if( length(fixConstr) )
    {
      if( parIndex %in% fixConstr[["indexes"]] )
      {
      constr <- lapply(fixConstr, function(input) input[parIndex])
      constr[["indexes"]]  <- 1
      } else {
        constr <- list()
      }
    }
    
    # Second call to .clusterSetUp (it's here for a reason!)
    if( multicore ){ 
      # Force evaluation of everything in the environment, so it will available to singleChain on cluster
      .forceEval(ALL = TRUE)
      
      tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    }
    
    grads[[parNam]] <- withCallingHandlers({
      laply(ranges[[parNam]],
            function(inCov, ...) { 
              
              sapply(1:nreps, 
                     function(input, ...) synGrad(object = object, 
                                                  param = fixParam[parNam], 
                                                  nsim = nsim, 
                                                  covariance = diag(inCov, 1)^2, 
                                                  addRegr = addRegr, 
                                                  nBoot = nBoot, 
                                                  constr = constr, 
                                                  ...)$gradient, 
                     ...)
              
            }, 
            .parallel = multicore,
            .inform = verbose, 
            ...)
    }, warning = function(w) {
      # There is a bug in plyr concerning a useless warning about "..."
      if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
    
    rownames(grads[[parNam]]) <- ranges[[parNam]]
    
  }
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  bestCov <- numeric( length(ranges) )
  names(bestCov) <- names(ranges)
  
  # Plotting signal the gradients for each parameter and variance
  par(mfrow = c(2, 1))
  for(parNam in names(ranges))
  {
    tmp <- t(grads[[parNam]])
    normGrad <- sweep(tmp, 2, apply(tmp, 2, sd), "/")
    boxplot( tmp, main = paste(parNam, ": gradient"), 
             ylab = paste("df/d", parNam, sep = ''))
    abline(h = 0, lty = 2)
    
    boxplot( sign(normGrad) * sqrt(abs(normGrad)), main = paste(parNam, ": scaled gradient"), 
             ylab = paste("(df/d", parNam,") / sdev(df/d", parNam, ")", sep = ''))
    abline(h = 0, lty = 2)
    
    bestCov[parNam] <- as.numeric( readline(prompt = paste("Enter the chosen standard deviation for ", 
                                                           parNam, 
                                                           " (probably in [", 
                                                           min(ranges[[parNam]]), ", ", max(ranges[[parNam]]),
                                                           "]) \n", sep = '')) )
  }
  
  return( invisible(list("bestCov" = bestCov, "grads" = grads)) )
}
