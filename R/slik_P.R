#############
#' Evaluates the synthetic log-likelihood.
#' 
#' @param param Vector of parameters at which the synthetic likelihood will be evaluated.
#' @param nsim  Number of simulation from the model.  
#' @param simulator function simulate data  
#' @param d desgin         
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
slik_P <- function(param, nsim, obsData, simulator,summaryStat,d,uRN) 
{
  .slik_P(param = param, 
        nsim  = nsim,
        obsData = obsData,
        simulator = simulator,
        summaryStat = summaryStat,
        d=d,
        uRN=uRN
        )$llk
  
}

####
# INTERNAL
####

.slik_P <- function(param, nsim, obsData, simulator,summaryStat,d,uRN) 
{
  
  if( !is.vector(param) ) stop("param should be a numeric vector.")
  
  simulData <- simulator(param,nsim,d,uRN) 
  
  simulSS <- summaryStat(simulData,d)
  
  simulSS <- .clean(X = simulSS, verbose = verbose)$cleanX
  
  simulSS <- simulSS[is.finite(rowSums(simulSS)),] #### MAHASEN remove infinit values. 
  
  
  # if(nrow(simulSS) < nsim / 3) warning(paste(nsim - nrow(simulSS), "out of", nsim, "statistics vectors", "contain NAs and will not be used"))
  # 
  if(nrow(simulSS) < nsim / 3) stop(paste(nsim - nrow(simulSS), "out of", nsim, "statistics vectors", "contain NAs and will not be used"))
  

  # Transforming the observation into stats
  obsStats <- summaryStat(obsData,d)
   
  # Dealing with possibly multiple datasets: each row could be a set of summary statistics
  if( !is.matrix(obsStats) ) obsStats <- matrix(obsStats, 1, length(obsStats))
  
  # Calculating log-likelihood
  out <- list("llk" = demvn(y = obsStats, X = simulSS, log = TRUE, verbose = FALSE), "mix" = 0)
 
  out$llk <- sum(out$llk)
  
  return( out )
  
}


