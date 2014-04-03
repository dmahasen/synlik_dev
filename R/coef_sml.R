#########
#' Parameters estimates for "sml" objects
#'
#' @param object  ("sml") object.
#' @param lag     (integer) final estimate of the parameter is the mean of the 
#'                last "lag" iterations. Hence lag < object@@niter.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @method coef sml
#' @export
#' 

coef.sml <- function(object, burn = 0, ...)
{
  if( burn > 0 ) object@estim <- object@estim[-(1:burn), ] 
  
  est <- colMeans( object@estim )
  
  names(est) <- names(object@param)
  
  return(est)
}

setMethod("coef",
          signature = signature(object = "sml"),
          definition = coef.sml
)