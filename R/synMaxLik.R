########
### Updates the parameters
########
## Does a Newton-Raphson update of the parameters
## maxJump is a scalar that is used to limit the size of the maximum move

.updateParams <- function(currPar, 
                          grad, 
                          covar,
                          gain, 
                          constr, 
                          quant)
{
  stopifnot( is.vector(currPar), is.vector(grad), is.matrix(covar) )
  
  nPar <- length( grad )
  
  # Propose a jump and limit its length so that it doesn't go
  # behond a certain quantile of the local normal model
  delta <- drop( gain * covar %*% grad )
  
  lim <- qchisq(quant, nPar)
  quad <- drop( crossprod(delta, qr.solve(covar, delta)) )
  delta <- delta * min( 1, sqrt(lim / quad) )
  
  newPar <- currPar - delta
  
  # Appling contraints to the new parameter vector
  if(length(constr))
  {
    newPar[constr$indexes] <- pmin(newPar[constr$indexes], constr$upper)
    newPar[constr$indexes] <- pmax(newPar[constr$indexes], constr$lower)
  }
  
  return( newPar )
  
}

# Estimated the hessian using a moving average, which discards negative
# definite hessians and outliers
# Arguments:
# - initHess: initial Hessian, useful at initialization (when length(hess) < lag).
# - hess: list of hessian from past iterations.
# - lag: number of past hessians considered.
# - zmax: z value used to discard the hessians the are too distant from the median
#         hessian.
# - inLag: maximum number of "initHess" to be concatenated to "hess" if the latter
#          is shorter than "lag". Useful when you want to have large "lag", but you
#          do not want to patch "initHess" lot of time at the beginning.
.getHessian <- function(initHess, hess, lag, zmax, inLag)
{
  
  # Remove hessians that are not positive definite
  eig <- sapply(hess, function(mat) min(eigen(mat, only.values = TRUE)$values)) 
  hess <- hess[ eig > 1e-14 ]
  
  hess <- tail(hess, min(lag, length(hess)))
  
  # Are we still in an early iteration (iter < lag)? 
  # If yes, we might have enough hessians, the first one is repeated npatch times
  npatch <- max(lag - length(hess), 0) 
    
  if( npatch ){ hess <- c(rlply(min(npatch, inLag), {initHess}), hess) }
    
  hess <- abind(hess, along = 3)
  
  # Median hessian (elementwise)
  med <- apply(hess, c(1, 2), median) 
  
  # Log sum of absolute of (elemetwise) differences from median Hessian 
  dist <- abs( sweep(hess, 1:2, med) )
  dist <- sqrt( apply(dist, 3, sum) )
  
  # Testing which hessian are outliers           
  dist_m <- median(dist)                   # Median
  dist_sd <- median( abs(dist - dist_m) )  # MAD
  
  # Calculate the indexes of the hessians to be used in the moving average
  # If(dist_sd == 0) the MAD is zero, hence the median hessian initial one
  # (1:length(z) <= npatch) is used to keep using the initial hessian as much
  # as we can, for robustness.
  good <- if(dist_sd == 0){ 
    
    which( dist == 0 )
    
  } else {
    
    z <- abs( (dist - dist_m) / dist_sd )
    
    which( z < zmax | (1:length(z) <= npatch) )
    
  }
  
  print(good)
   
  out <- apply(hess[ , , good], 1:2, mean)
  
  return( out )
  
}

# Test for .getHessian
# library(plyr)
# library(MCMCpack)
# 
# m <- rlply(100, rwish(100, matrix(c(1,.3,.3,1), 2, 2)))
# 
# m[[101]] <- rwish(100, matrix(c(2,.7,.7,1.2), 2, 2))
# 
# .getHessian(hess = m, lag = 10, zmax = 3)
# 
# m <- rlply(5, rwish(100, matrix(c(1,.3,.3,1), 2, 2)))
# 
# m[[6]] <- rwish(100, matrix(c(2,.7,.7,1.2), 2, 2))
# 
# .getHessian(hess = m, lag = 10, zmax = 3)



#######
#### Gets positive-definite covariace out of hessian 
#######
#
# ARGS:
# hessian = (matrix)
# uplim and uplim = (numeric) of length ncol(hessian) that limit the size of the diagonal
#                   elements of the covariance                    
# verbose = (logical) print whether the limits are violated

.getCovariance <- function(hessian, upLim, lowLim, verbose = FALSE) 
{
  stopifnot( is.matrix(hessian), is.vector(upLim), is.vector(lowLim), all(upLim > lowLim))
  
  nPar <- nrow(as.matrix(hessian))
  
  # Inverting Hessian to get the covariance
  covar <- .qrInverse(hessian, tolQR = 0, imposePD = TRUE, tilt = 1e-5)
  
  if( any(diag(covar) > upLim) || any(diag(covar) < lowLim) )
  {
    if(verbose == TRUE){ 
      print("Variance of proposal is too high or too low! (probably to high)") 
      print(diag(covar))
    }
    
    corr <- cov2cor( covar ); 
    sdev <- sqrt(diag(covar))
    sdev <- pmax( pmin(sdev, sqrt(upLim)), sqrt(lowLim))
    covar  <- diag(sdev, nPar) %*% corr %*% diag(sdev, nPar)
  }
  
  return( covar )
}

#########
#### The stochastic optimization routine
#########
#' Synthetic likelihood maximization routine.
#'
#' @param object  ("synlik") object.
#' @param nIter   (integer) numer of iterations.
#' @param nsim    (integer) numer of simulations from the model at each step.
#' @param initCov  (matrix) initial covariance matrix used to simulate the paramters at each step.
#' @param initPar  (numeric) vector of initial values of the parameters.
#' @param addRegr  (logical) 
#'           if FALSE the statistics calculated by object@@summaries will be used (SL approach)
#'           if TRUE the simulated parameters will be regressed on the statistics and the 
#'           fitted values of the paramaters given the _observed_ statistics will be used as statistics
#'            (SL+ approach)
#' @param constr (named list) of 3 elements:
#'           [["indexes"]] = (numeric integers) indexes of the elements to check;
#'           [["upper"]]  = (numeric) upper bounds for the elements in "indexes";
#'           [["lower"]]  = (numeric) lower bounds for the elements in "indexes".
#' @param control Named list of control setting for the optimization routine.
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param verbose  (logical) if TRUE lots of things will be printed.
#' @param ...  additional arguments to be passed to object@@simulator and object@@summaries.
#'             In general I would avoid using it and including in those two function everything they need.
#' @return object of class "synOptim": see info about that class.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @export
 
synMaxlik <- function(object, nIter, nsim, 
                      initCov, initPar = object@param, 
                      nBoot = 0,
                      adapt = FALSE,
                      fixVar = TRUE,
                      addRegr = TRUE, constr = list(), control = list(),
                      multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                      verbose = FALSE,  ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  stopifnot(is.vector(initPar), is.matrix(initCov))
  
  nPar <- length(initPar)
  currCov <- oldCov <- initCov
  currPar <- oldPar <- initPar
  
  # Setting up control parameter
  ctrl <- list( "gain" = 1, 
                "gainExp" = 0.602, 
                "lag" = max(round(nIter / 3), 1),
                "limCov" = list("upper" = diag(currCov) * Inf, 
                                "lower" = diag(currCov) * 0),
                "quant" = 0.25,
                "zmax" = 2, 
                "inLag" = max(round(nIter / 3), 1),
                "initHess" = solve(currCov),
                "recycle" = FALSE)
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if(!identical(names(ctrl$limCov), c("upper", "lower"))) stop("limCov should have names \"upper\" and \"lower\" (in that order)")
  
  # Checking if the contraints are valid
  if( length(constr) )
    stopifnot( all(c("indexes", "upper", "lower") %in% names(constr)),
               length(constr$upper) == length(constr$lower), length(constr$lower) <= nPar, 
               length(constr$indexes) == length(constr$lower), all(constr$upper > constr$lower)
               ) 
    
  # Create decreasing gain sequence
  gainSeq <- ctrl$gain / ( (1:nIter) ^ ctrl$gainExp )
  
  # I will store the results here
  resultPar <- resultGrad <- matrix(NA, nIter, nPar)
  resultHess <- resultCovar <- list()
  resultLoglik <- numeric(nIter) 
  
  # Setting up a cluster if needed
  if(multicore) {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik")
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
  }
  
  # List where recycled values, parameters and weights will be stored
  storage <- list()
  
  # The main loop of the optimization
  for(ii in 1:nIter)
  {
    # Calculate gradient and hessian                     
    tmp <- synGrad(object, param = currPar, nsim = nsim, covariance = currCov, 
                   addRegr = addRegr, fixVar = fixVar, nBoot = nBoot, constr = constr, 
                   multicore = multicore, ncores = ncores, cluster = cluster, ...)  
      
    gradient <- resultGrad[ii, ] <- - tmp$gradient
    resultHess[[paste("Iter", ii, sep = "")]] <- - tmp$hessian
    resultLoglik[ii] <- tmp$llk
    storage[[ii]] <- tmp[["stored"]]
    
    # Get Hessian from moving average
    hessian_hat <- .getHessian(ctrl$initHess, resultHess, lag = ctrl$lag, zmax = ctrl$zmax, inLag = ctrl$inLag)
    
    # Get the covariance which will be used to simulate the parameters at the next iteration
    if( adapt )
    {
      currCov <- .getCovariance(hessian = hessian_hat, upLim = ctrl$limCov$upper, lowLim = ctrl$limCov$lower)
    } 
   
    # Update parameters values
    currPar <- .updateParams(currPar = currPar, grad = gradient, 
                             covar = currCov, 
                             gain = gainSeq[ii], 
                             constr = constr, quant = ctrl$quant)
    
    oldCov <- currCov
    oldPar <- currPar
        
    resultPar[ii, ] <- currPar
    resultCovar[[paste("Iter", ii, sep = "")]] <- currCov
    
    if(verbose){ 
      print( paste("Iteration", ii, ", parameters =", currPar ) )
      print( currCov )
    }
    
  }
  
  if(multicore && clusterCreated) stopCluster(cluster)
  
  colnames(resultPar) <- names(object@param)
  for(ii in 1:nIter) dimnames(resultHess[[ii]]) <- dimnames(resultCovar[[ii]]) <- list(names(object@param),
                                                                                       names(object@param))
  
  # Setting up control list for "continue.synMaxlik" method
  toContinue <- ctrl
  toContinue$initHess <- hessian_hat
  toContinue$gain <- tail(gainSeq, 1) 
  
  ### Class Definition
  .synMaxlik( object,
              initPar = initPar,
              niter   = as.integer(nIter),
              nsim    = as.integer(nsim),
              initCov = initCov,
              addRegr  = addRegr,
              constr = constr,
              control = ctrl,
              continueCtrl = toContinue,
              multicore = multicore,
              ncores = as.integer(ncores),
                                        
              resultPar = resultPar,
              resultGrad = resultGrad,
              resultHess = resultHess,
              resultCovar = resultCovar,
              resultLoglik = resultLoglik, 
              storage = storage)
  
}


