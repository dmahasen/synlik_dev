#######
### Method to estimate the parameters through MCMC
#######
#' MCMC parameter estimation for objects of class \code{synlik}.
#' 
#' @param object An object of class \code{synlik}.
#' @param initPar see \code{\link{smcmc-class}}.
#' @param niter see \code{\link{smcmc-class}}.
#' @param nsim  see \code{\link{smcmc-class}}.
#' @param propCov see \code{\link{smcmc-class}}.
#' @param burn see \code{\link{smcmc-class}}.
#' @param nchains number of MCMC chains used, by default 1.
#' @param priorFun see \code{\link{smcmc-class}}.
#' @param targetRate see \code{\link{smcmc-class}}.
#' @param recompute see \code{\link{smcmc-class}}.             
#' @param multicore  see \code{\link{smcmc-class}}.
#' @param ncores   see \code{\link{smcmc-class}}.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param control see \code{\link{smcmc-class}}.
#' @param ... additional arguments to be passed to \code{slik} function, see \code{\link{slik}}.
#' @return If \code{nchains == 1} an object of class \code{smcmc} will be returned. Otherwise a list of \code{nchains}
#'         \code{smcmc} objects will be returned.
#'
#' @details When \code{multicore == TRUE} there are two possible behaviours: if \code{nchains == 1} the 
#'          the parallelism takes place at each step of the MCMC chain (\code{object@@simulator} and 
#'          \code{object@@summaries} are called in parallel). If \code{nchains > 1} each chain will be assigned
#'          to a different node of the cluster, and there is no parallelism within the chains. 
#'         
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>, code for adaptive step from the adaptMCMC package.   
#' @references Vihola, M. (2011) Robust adaptive Metropolis algorithm with coerced acceptance rate. 
#'             Statistics and Computing.           
#' @export
#' 
smcmc <- function(object, 
                  initPar, 
                  niter, 
                  nsim,
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
  
  # Force evaluation of everything in the environment, so it will available to funToApply on cluster
  if( multicore ) .forceEval(ALL = TRUE)
  
  # multicore, ncores and cluster go in the ...
  funToApply <- function(notUsed, ...)
  {
    .smcmc(object = object, 
           initPar = initPar, 
           niter = niter, 
           nsim = nsim,
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
    if(multicore)
    {
      tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
      cluster <- tmp$cluster
      ncores <- tmp$ncores
      clusterCreated <- tmp$clusterCreated
    }
        
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
  
  if(nchains == 1) out <- out[[1]]
  
  return(out)
}
