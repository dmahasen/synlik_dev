
#####
#' Cross-validation for the empirical saddlepoint density
#' @description Performs cross-validation to choose the \code{decay} tuning parameters 
#'              which determines the mixture between a consistent and a Gaussian estimator
#'              of the Cumulant Generating Function (ECGF).
#'
#' @param decay Numeric vector containing the possible value of the tuning parameter.
#' @param simulate  Function with prototype \code{function(...)} that will be called \code{nrep} times
#'                  to simulate \code{d}-dimensional random variables. 
#'                  Each time \code{simulator} is called it will return a \code{n} by \code{d} matrix.
#' @param nrep Number of times the whole cross-validation procedure will be repeated, by calling
#'             \code{simulator} to generate random variable and computing the cross-validation score
#'             for every element of \code{decay}.
#' @param multicore  see \code{\link{smcmc-class}}.
#' @param ncores   see \code{\link{smcmc-class}}.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param control A list of control parameters, with entries:
#'         \itemize{
#'         \item{ \code{K} }{The number of folds to be used in cross-validation. By defaults \code{K = 10};}
#'         \item{ \code{mixMethod} }{ Method used to estimate the ECGF. By default equal to "mse". See \code{\link{ecgf}} for details; }
#'         \item{ \code{draw} }{ If \code{TRUE} the results of cross-validation will be plotted. \code{TRUE} by default;}
#'         \item{ \code{nNorm} }{ Number of simulations to be used in order to estimate the normalizing constant of the saddlepoint density.
#'                                By default equal to 1e3.}
#'         }
#' @return A list with entries:
#'         \itemize{
#'          \item{ \code{summary} }{ A matrix of summary results from the cross-validation procedure.  }
#'         \item{ \code{negLogLik} }{A matrix \code{length{decay}} by \code{control$K*nrep} where the i-th row represent the negative loglikelihood
#'                                   estimated for the i-th value of \code{decay}, while each column represents a different fold and repetition.}
#'         \item{ \code{normConst} }{ A matrix \code{length{decay}} by \code{nrep} where the i-th row contains the estimates of the normalizing constant.}
#'         }
#'         is returned invisibly.
#'         If \code{control$draw == TRUE} the function will also plot the cross-validation curve. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' 
#' # Saddlepoint is needed
#' selectDecay(decay = seq(0.05, 0.1, 0.5, 1), 
#'             simulator = function(...) rgamma(100, 2, 1), 
#'             nrep = 4)
#'             
#' # Saddlepoint is not needed
#' selectDecay(decay = seq(0.05, 0.1, 0.5, 1), 
#'             simulator = function(...) rnorm(100, 0, 1), 
#'             nrep = 4)
#' 
#' @export
#'

selectDecay <- function(decay, 
                        simulator, 
                        nrep, 
                        multicore = !is.null(cluster),
                        cluster = NULL,
                        ncores = detectCores() - 1, 
                        control = list(), 
                        ...)
{
  
  # Control list which will be used internally
  ctrl <- list( "K" = 10,
                "mixMethod" = "mse",
                "draw" = TRUE,
                "nNorm" = 1000 )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available on cluster
    .forceEval(ALL = TRUE)
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  #  We simulate data
  datasets <- rlply(nrep, simulator(...))
    
  # Get cross-validated negative log-likelihood for each dataset
  withCallingHandlers({
    tmp <- llply(datasets, .selectDecay, 
                 .progress = "text", .parallel = multicore, 
                 # Extra args for .selectDecay
                 decay = decay, K = ctrl$K, mixMethod = ctrl$mixMethod, nNorm = ctrl$nNorm) 
  }, warning = function(w) {
    # There is a bug in plyr concerning a useless warning about "..."
    if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
      invokeRestart("muffleWarning")
  })
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  negLogLik <- do.call("cbind", lapply(tmp, "[[", "negLogLik"))
    
  # Preparing output: returning a summary table of scores, all -log-lik and all normalizing constants
  out <- list()
  
  out$negLogLik <- do.call("cbind", lapply(tmp, "[[", "negLogLik"))
  colnames(out$negLogLik) <-  1:(ctrl$K*nrep)
  
  out$summary <- rbind( decay, rowMeans(out$negLogLik), apply(out$negLogLik, 1, sd) )  
  rownames(out$summary) <- c("decay", "mean_score", "sd_score")
  
  out$normConst <- do.call("cbind", lapply(tmp, "[[", "normConst"))
  rownames(out$normConst) <- decay
  
  # (Optionally) Plot the score for each value in "decay".
  if(ctrl$draw){
    par(mfrow = c(2, 1))
    
    .forestPlot( out$summary[2, ], out$summary[3, ], xlabels = decay, xlab = "Decay", ylab = "- Log-likelihood (mean+-sd)" )
    
    .forestPlot( rowMeans(out$normConst), apply(out$normConst, 1, sd), xlabels = decay, 
                 ylab = "Normalizing constants (mean+-sd)", xlab = "Decay" )
    
  }
    
  return( invisible( out ) )
  
}



########
# Internal R function
########

.selectDecay <- function(X, decay, K, mixMethod, nNorm)
{
  if( !is.matrix(X) ) X <- matrix(X, length(X), 1)
  
  n <- nrow(X)
  
  # Divide the sample in K folds
  folds <- mapply(function(a, b) rep(a, each = b), 1:K, c(rep(floor(n / K), K - 1), floor(n / K) + n %% K), SIMPLIFY = FALSE )  
  folds <- do.call("c", folds)
  
  ngrid <- length(decay)
  
  negLogLik <- matrix(NA, ngrid, K)
  normConst <- numeric(ngrid)
  rownames(negLogLik) <- decay
  colnames(negLogLik) <- 1:K
  
  sam <- rmvn(nNorm, colMeans(X), cov(X))
  # For each value of "decay" and for each fold, calculate the negative log-likelihood 
  # of the sample points belonging to that fold.
  for(ii in 1:ngrid)
  {
    
    normConst[ ii ] <- mean( dsaddle(y = sam, X = X, decay = decay[ii], mixMethod = mixMethod, log = FALSE)$llk /  dmvn(sam, colMeans(X), cov(X)) )
    
    negLogLik[ii, ] <- aaply(1:K, 
                             1,
                             function(input){
                               index <- which(folds == input)
                               -sum( dsaddle(X[index, , drop = F], X = X[-index, , drop = F], 
                                             decay = decay[ii], mixMethod = mixMethod, log = TRUE)$llk ) + length(index) * log(normConst[ ii ])
                             })
    
  }
  
  return( list( "negLogLik" = negLogLik, "normConst" = normConst ) )
  
}
