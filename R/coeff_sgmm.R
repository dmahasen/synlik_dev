

coeff.sgmm <- function(object, nBurn, 
                       fixVar = FALSE,
                       multicore = !is.null(cluster),
                       cluster = NULL,
                       ncores = detectCores() - 1, 
                       verbose = FALSE,
                       ...)
{
  # Extract simulated parameters (param, X) and statistics (y), 
  # after discarding those simulated in the first nBurn steps.
  param <- lapply(object@storage[-(1:nBurn)], "[[", "X")
  param <- do.call("rbind", param)
  
  y <- lapply(object@storage[-(1:nBurn)], "[[", "y")
  y <- do.call("rbind", y)
  
  # Transform observations into summary statistics
  summaries <- object@summaries
  yobs <- if( !is.null(summaries) ) summaries(x = object@data, extraArgs = object@extraArgs, ...) else object@data
  
  # Call function that does the smooth and give you back the objective function
  objFun <- smoothLik(yobs = yobs, y = y, param = param, fixVar = fixVar, multicore = multicore,
                      cluster = cluster, ncores = ncores, verbose = verbose)
  
  return( objFun )
}

# TESTED

# yobs = observed statistics (vector)
# y    = simulated statistics (matrix)
# param = simulated parameters (matrix)
smoothLik <- function(yobs, y, param,
                      fixVar = FALSE,
                      multicore = !is.null(cluster),
                      cluster = NULL,
                      ncores = detectCores() - 1, 
                      verbose = FALSE, 
                      ...)
{
  if( !is.matrix(y) ) y <- matrix(y, length(y), 1)
  
  np <- ncol(param)
  n <- nrow(param)
  ns <- ncol(y)
  
  X <- param
  
  # Construct regression matrix for mean model (intercept + linear + quadratic + interactions) 
  Xm <- .quadModMat(X)
  colnames(Xm) <- sapply(1:ncol(Xm), function(input) paste("x", input, sep = ''))
  Xm <- as.data.frame(Xm)
  
  ### Construct regression matrix for variance model (no intercept!! + linear) 
  Xv <- scale(X)
  center <- attr(Xv, "scaled:center")
  scale <- attr(Xv, "scaled:scale")
  colnames(Xv) <- sapply(1:ncol(Xv), function(input) paste("x", input, sep = ''))
  Xv <- as.data.frame(Xv)
  
  weights <- rlply(ns, NULL) 
  
  funToApply <- function(ii, ...) {
    # Building regression formula: 1 + s(par_1) + ... + s(par_np) + interactions
    form <- list()
    form[[1]] <- paste("y[ , ", ii, "] ~ 1", sep = '')
    form <- append(form, lapply(2:(np+1), function(input) paste("+ s(x", input, ")", sep = '')))
    if(np > 1) form <- append(form, lapply((2+2*np):ncol(Xm), function(input) paste("+ x", input, sep = '')))
    form <- as.formula( do.call("paste", form) )
        
    return( gam(form, data = Xm, weights = weights[[ii]], ...) )
  }
  
  if( multicore ){ 
    verbose <- FALSE
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = c("synlik", "mgcv"))
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
    
    environment(funToApply) <- .GlobalEnv
  }
    
  # In the first loop the smooth regressions for the mean of each statistic 
  # assume homoschedasticity, in the second loop I use weights equal to the
  # reciprocal of the root estimated variance.
  for(kk in 1:2)
  { 
    
    if(multicore) { clusterExport(cluster,  list("y", "Xm", "weights", "np", "ns"), envir = environment()) }
    
    # Loop over the statistics to estimate the mean of each
    withCallingHandlers({
      regr <- llply(1:ns, 
                    funToApply,
                    .parallel = multicore,
                    .inform = verbose, 
                    ...)
    }, warning = function(w) {
      # There is a bug in plyr concerning a useless warning about "..."
      if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
    
    # Extract residuals from each regression, square and log them.
    res <- t( laply(regr, "[[", "residuals") )
    if(ncol(res) > nrow(res)) res <- t(res) 
    sqRes <- res ^ 2
    
    if( !fixVar )
    {
    # Regress squared-residuals on parameters
    fixCov <- robCov( t(res) )
    varCoef <- expRegr(y = sqRes, X = Xv, weights = fixCov$weights)
    
    # Calculate weights for 2nd mean regression: 1 / sqrt( var_hat )
    fitVar <- exp( cbind(1, as.matrix(Xv)) %*% varCoef )
    weights <- alply(fitVar,
                     2,
                     function(input) {  
                       w <- 1 / sqrt(input)
                       w <- w * (n / sum(w))
                       } )
    } else {
      break;
    }
  }
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  # Variance of the residuals will be needed by objFun (which is a closure)
  fixCov <- robCov( t(res) )
  fixDec <- robCov(t(res) / apply(res, 2, sd)) 
  
  # Synthetic Likelihood function based on the estimated model for the mean and 
  # the marginal variances of the statistics wrt the parameters.
  objFun <- function(param, fun = "sl")
  {
    nX <- drop( param )
    
    # Transform the input parameters into a data frame to use predict.gam()
    nXm <- .quadModMat( matrix(nX, 1, length(nX)) )
    colnames(nXm) <- names(Xm)
    nXm <- as.data.frame(nXm)
    
    # Mean (mu) of the statistics at "nX", predicted by gam()
    mu <- sapply(regr, predict, newdata = nXm)
    
    # Marginal standard deviations (sdev) of the statistics at "nX", predicted by expRegr() or constant
    if( !fixVar )
    {
    nXv <- cbind(1, scale(t(nX), center, scale))
    sdev <- drop( exp(0.5 * nXv %*% varCoef) )
    } else {
      sdev <- fixCov$sd
    }
            
    # Return either synthetic likelihood or gmm function
    out <- sum( (fixDec$E %*% diag(sdev^-1, ns) %*% drop(yobs - mu))^2 )
    
    if(fun == "sl") out <- - out/2 - fixDec$half.ldet.V - sum(log(sdev)) - log(2 * pi) * length(yobs) / 2
        
    return( out )
    
  }
  
  return( objFun )
  
}


