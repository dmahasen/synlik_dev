########
#' Log-likelihood of \code{smcmc} objects
#'
#' @param object  ("smcmc") object.
#' @param nreps   (integer) number of simulations used to tilt the negative Hessian (-H)
#'                toward positive definiteness (PD). Used only if the -H is not PD and \code{method == "asym"}.
#' @param boot    (logical) relevant only if -H is not PD. 
#'                If TRUE hessians will be simulated by resampling parameters and likelihoods. 
#'                If FALSE hessians will be simulated the asymptotic distribution of the regression
#'                coefficients.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @method logLik smcmc
#' @export
#' 

logLik.smcmc <- function(object,
                         
                         burn = 0,
                         quant = 0.1,
                         
                         nreps = 1000, 
                         boot = TRUE, 
                         ...)
{ 
  # Including only parameters that were estimated
  varPar <- diag(object@propCov) > 0
  
  mcmcObj <- list("initPar" = object@initPar,
                  "propCov" = object@propCov, 
                  "llkStore" = object@llkStore,
                  "parStore" = object@parStore)
  
  llkMax <- logLik.mcmc(object = mcmcObj,
                        nreps = nreps, 
                        boot = boot, 
                        quant = quant,
                        burn = burn,
                        ...)
  
  return(llkMax) 
}

setMethod("logLik",
          signature = signature(object = "smcmc"),
          definition = logLik.smcmc
)