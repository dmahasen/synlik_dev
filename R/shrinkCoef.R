##########################
#' Get coefficients of ridge regression or lasso
#' 
#' @description This function simulates parameters in an hypercube and uses these parameters sets
#'              to simulated statistics from the model. The simulated values of each parameter 
#'              are then ridge-regressed on the simulated statistic. The function finds the optimal
#'              penalty and return a set of regression coefficients.
#' 
#' @param object An object of class "synlik".
#' @param nsim Number of summary statistics to be simulated.
#' @param mu   mean around which the parameters will be simulated.
#' @param sigma covariance matrix used to simulate the parameters.
#' @param ... additional arguments to be passed to object@@simulator and object@@summaries.
#'            In general I would avoid using it and including in object@@extraArgs everything they need.
#' @return A list where [["ridgeCoef"]] is a matrix where the i-th row contains the ridge-regression coefficient
#'         resulting from regressing the i-th parameter on the statistics; [["meanStats"]] is a vector containing
#'         the means of the simulated statistics; [["sdevStats"]] is a vector containing the standard deviations of the 
#'         simulated statistics; [["penalty"]] is a vector of optimal penalties from each ridge regression.
#'            
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>, re-using much code from \code{MASS::lm.ridge()}.                         
#' @export shrinkCoef       

shrinkCoef <- function(object, nsim, mu, sigma, type = "ridge", constr = list(), 
                       multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                       verbose = TRUE, clean = TRUE, ...)
{
  if( !(type %in% c("ridge", "lasso")) ) stop(" type should be either \"ridge\" or \"lasso\" ")
  
  nPar <- length(object@param)
  
  stopifnot( is.numeric(mu), is.numeric(sigma), length(mu) == nPar )
  
  # Fixing parameters with variance == 0
  fixPar <- which( diag(sigma) == 0 )
  if( length(fixPar) ) message( paste("Parameters", names(object@param)[fixPar], "are considered fixed") )
  
  # Simulate the parameters from truncated multivariate normal
  param <- .paramsSimulator(theMean = mu, covar = sigma, nsim = nsim, constr = constr)
  
  if( multicore ){ 
    # Force evaluation of everything in the environment
    .forceEval(ALL = TRUE)
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  # Simulate a vector of statistics for each parameter value
  withCallingHandlers({
    simul <- alply(.data = param,
                   .margins = 1,
                   .fun = function(input){
                     
                     return(
                       tryCatch( .simulate.synlik(object = object, nsim = 1, param = input, stats = TRUE, clean = FALSE, multicore = FALSE, ...),
                                 error = function(e) return(NA))
                     )
                     
                   },
                   .parallel = multicore,
                   ...
    )
  }, warning = function(w) {
    # There is a bug in plyr concerning a useless warning about "..."
    if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
      invokeRestart("muffleWarning")
  })
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  simul <- do.call("rbind", simul)
  
  # Removing missing values
  if(clean)
  {
    goodPar <- !is.na(simul[ , 1])
    simul <- simul[goodPar, , drop = FALSE]
    param <- param[goodPar, , drop = FALSE]  
  }
  
  regrCoef <- matrix(NA, nPar, ncol(simul) + 1) # + 1 is intercept
  penalty <- numeric(nPar)
  names(penalty) <- rownames(regrCoef) <- names(object@param)
  
  # For each parameter, find the optimal penalty and parameters using ridge regression or lasso
  for(iPar in 1:nPar)
  {
    if( !(iPar %in%fixPar) )
    {
      y <- param[ , iPar]
      
      if(type == "ridge")
      {
        tmp <- .autoLm.ridge(y ~ simul)
        
        # Put the coefficients on the original scale
        if( is.null(tmp$Inter) ) stop("For some reason the ridge regression doesn't have an intercept!")
        
        scaledCoef <- t(as.matrix(tmp$coef / tmp$scales))
        inter <- tmp$ym - scaledCoef %*% tmp$xm
        
        regrCoef[iPar, ] <- c(inter, scaledCoef)
        penalty[iPar] <- tmp$lambda
      } 
      
      if(type == "lasso")
      {
        tmp <- .autoLm.lasso(x = simul, y = y, verbose = verbose)
        
        if(verbose) title( main = names(object@param)[iPar] )
        
        # Store the regression coefficients with the optimal penalty
        regrCoef[iPar, ] <- tmp$coef 
        penalty[iPar] <- tmp$fraction 
      }
      
    }
    
  }
  
  return( list("regrCoef" = regrCoef, "penalty" = penalty) )
}
