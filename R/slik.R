#############
#' Evaluates the synthetic log-likelihood.
#' 
#' @param object An object of class \code{synlik}.
#' @param param Vector of parameters at which the synthetic likelihood will be evaluated.
#' @param nsim  Number of simulation from the model.             
#' @param multicore  (logical) if \code{TRUE} the \code{object@@simulator} and \code{object@@summaries} functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if \code{multicore == TRUE}.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to \code{object@@simulator} and \code{object@@summaries}.
#'            In general I would avoid using it and including \code{object@@extraArgs} everything they need.
#' @return The estimated value of the synthetic log-likelihood at \code{param}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>    
#' @references Simon N Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102--1104, 2010.
#' @examples
#' data(ricker_sl)
#' set.seed(643)
#' slik(ricker_sl, param = c(3.8, -1.2, 2.3), nsim = 500)                     
#' @export
#' 
slik <- function(object, param, nsim, multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, ...) 
{
 
  .slik(object = object, 
        param = param, 
        nsim  = nsim, 
        decay = 0.5, 
        multicore = multicore, 
        ncores = ncores, 
        cluster = cluster, 
        ...)$logLik
  
}

####
# INTERNAL
####

.slik <- function(object, param, nsim, saddle = FALSE, decay = 0.5, multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, ...) 
{
  
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  if( !is.vector(param) ) stop("param should be a numeric vector.")
  
  # Simulating from model
  if(multicore)
  {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    
    coresSchedule <- c( rep(floor(nsim / ncores), ncores - 1), floor(nsim / ncores) + nsim %% ncores)
    
    simulData <- clusterApply(
      cluster, 
      coresSchedule, 
      function(input, ...)
      {
        .simulate.synlik(object, param = param, nsim = input, stats = TRUE, clean = TRUE, verbose = FALSE, ...)
      } 
      , ...
    )
    
    simulData <- do.call(rbind, simulData)
    
    if(clusterCreated) cluster <- stopCluster(cluster)
    
  } else {
    simulData <- .simulate.synlik(object, param = param, nsim = nsim, stats = TRUE, clean = TRUE, verbose = FALSE, ...) 
  }
  
  if(nrow(simulData) < nsim / 3) warning(paste(nsim - nrow(simulData), "out of", nsim, "statistics vectors", "contain NAs and will not be used"))
  
  # Transforming the observation into stats
  summaries <- object@summaries
  obsStats <- if( !is.null(summaries) ) summaries(x = object@data, extraArgs = object@extraArgs, ...) else object@data
  
  # Calculating log-likelihood
  return( .empDens(y = obsStats, X = simulData, saddle = saddle, decay = decay, tol = 1e-6, log = TRUE) )
  
}


