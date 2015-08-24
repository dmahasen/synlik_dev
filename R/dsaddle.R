#####
#' Empirical saddlepoint density 
#' @description Gives a pointwise evaluation of the empirical saddlepoint and optionally of
#'              its gradient at position y
#'
#' @param y Point at which the SPA is evaluated (d dimensional vector).
#' @param X n by d matrix containing the data.
#' @param tol Tolerance used to assess the convergence of the rootfinding routine used to fit
#'            the saddlepoint density. Default value is 1e-6.
#' @param decay Rate at which the SPA falls back on a normal density. Should be a positive number,
#'              by default set to 0.5.
#' @param deriv If TRUE also the gradient of the log-saddlepoint density is returned.
#' @param log If TRUE the log of the saddlepoint density is returned.
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{llk} }{The value of the empirical saddlepoint at y;}
#'         \item{ \code{mix} }{The mixture of normal-saddlepoint used (1 means only saddlepoint);}
#'         \item{ \code{grad} }{The gradient of the log-density at y (optional);}
#'         }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#' @export
#'

dsaddle <- function(y, X, tol=1e-6, decay = 0.5, deriv = FALSE, mixMethod = "mse", log = FALSE, 
                    multicore = !is.null(cluster), ncores = detectCores() - 1, cluster = NULL) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  if( !is.matrix(X) ) X <- matrix(X, length(X), 1)
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available on cluster
    .forceEval(ALL = TRUE)
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  d <- ncol(X)
  
  # Offsetting dimensionality, so decay stays pretty much at the same level for any d.
  decay <- decay / ( d ^ 2 )
  
  if( !is.matrix(y) && d > 1 ) y <- matrix(y, 1, d)
  
  # Pre-calculating covariance
  #preCov <- .robCov(t(X), alpha2 = 4, beta2 = 1.25)
  preCov <- .robCov(t(X), alpha = 10, alpha2 = 10, beta2 = 1.25)
  
  # If there are some statistics with zero variance we remove them
  if( length(preCov$lowVar) )
  {
    stop("The columns of X indexed", preCov$lowVar, "have zero variance.")
    #y <- y[-preCov$lowVar]
    #X <- X[ , -preCov$lowVar, drop = FALSE]
  }
  
  # Weighting the statistics in order to downweight outliers
  X <- preCov$weights * X
  
  # Divide saddlepoint evaluations between cores
  withCallingHandlers({
    tmp <- alply(y, 1, .dsaddle, .parallel = multicore,
                 # Args for .dsaddle()
                 X = X,
                 preCov = preCov,
                 tol = tol,
                 decay = decay,
                 deriv = deriv,
                 mixMethod = mixMethod,
                 log = log)}, warning = function(w) {
                   # There is a bug in plyr concerning a useless warning about "..."
                   if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
                     invokeRestart("muffleWarning")
                 })
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  out <- list( "llk" = sapply(tmp, "[[", "llk"),
               "mix" = sapply(tmp, "[[", "mix"),
               "grad" = if(deriv) t( sapply(tmp, "[[", "grad") ) else NULL
  )
  
  return( out )
  
}


##########
# INTERNAL
##########

.dsaddle <- cmpfun(function(y, X, preCov, tol=1e-6, decay = 0.5, 
                            deriv = FALSE, mixMethod = "mse", 
                            maxit = 100, log = FALSE) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  
  if( !is.vector(y) ) y <- as.vector(y)
  
  m <- length(y)
  
  if(!is.matrix(X)){
    if(m > 1){ 
      stop("Error: simulated data must be entered in matrix form")
    }else{
      X <- matrix(X, length(X), 1)
    }
  }
  
  n <- nrow(X)
  
  # Initial guess of the root is the solution to the Gaussian case
  # the gain is one step less of Newton on average.
  lambda <- drop( crossprod(preCov$E, preCov$E %*% (y - preCov$mY)) )
  
  rss <-  drop( sum( (preCov$E %*% (y - preCov$mY))^2 ) )
  
  # Choosing method to determine the mixture between gaussian and normal ecgf
  switch(mixMethod,
         "gaus" = {
           mix <- exp( - decay * rss / (sqrt(m) * 2) )
         },
         "mse" = {
           sadMse <- exp(rss) 
           gausMse <- rss + 0.5*rss^2 + 1
           mix <- (gausMse / sadMse) ^ decay
         },
         stop("\"mixMethod\" must be either \"gaus\" or \"mse\"") 
  )
  
  if(decay == Inf) mix <- 0.0
  
  b <- .ecgf(lambda, X, kum1 = preCov$mY, kum2 = preCov$COV, mix = mix, grad = 2, mixMethod = mixMethod)
  
  ## Newton loop to minimize K(lambda) - t(lambda)%*%y or solve dK(lambda) = y wrt lambda
  kk <- jj <- 1
  # Convergence test: see [con_test] below.
  while( any( abs(b$dK-y) / sqrt(diag(preCov$COV)) > tol ) && kk < maxit ) 
  { 
    kk <- kk + 1
    
    D <- diag(diag(b$d2K)^-0.5, nrow = m, ncol = m)
    
    # Try solve scaled linear system fast, if that doesn't work use QR decomposition.
    d.lambda <- - drop( D %*% tryCatch(solve(D%*%b$d2K%*%D, D%*%(b$dK-y), tol = 0), error = function(e) qr.solve(D%*%b$d2K%*%D, D%*%(b$dK-y), tol = 0)) )
    
    lambda1 <- lambda + d.lambda ## trial lambda
    
    b1 <- .ecgf(lambda1, X, kum1 = preCov$mY, kum2 = preCov$COV, mix = mix, grad = 2, mixMethod = mixMethod)
    if ( sum( abs(b1$d2K) ) == 0 ) return(NA) ## spa breakdown (possibly too tight)
    
    jj <- 1
    c1 <- 10^-4
    alpha <- 1
    rho <- 0.5
    ## Line search checking Arminjo condition and step halving
    while( ( b1$K - crossprod(lambda1, y) ) > 
             ( b$K - crossprod(lambda, y) ) + c1 * alpha * crossprod(d.lambda, drop(b$dK-y)) && jj < 50)  
    {
      jj <- jj + 1
      alpha <- alpha * rho
      d.lambda <- d.lambda * alpha
      lambda1 <- lambda + d.lambda
      b1 <- .ecgf(lambda1, X, kum1 = preCov$mY, kum2 = preCov$COV, mix = mix, grad = 2, mixMethod = mixMethod)
    }
    
    ## now update lambda, K, dK, d2K
    lambda <- lambda1 
    b <- b1
    
  } ## end of Newton loop
  ## lambda is the SPA lambda...
    
  if(kk > 50 || jj > 20) warning(paste("The convergence of the saddlepoint root-finding is quite slow! \n",
                                       "Outer root-finding Newton-Raphson:", kk, "iter \n",
                                       "Inner line search for Arminjo condition:", jj, "iter"))
  
  spa <- b$K - crossprod(lambda, y) - 0.5 * log( 2*pi ) * m - 0.5 * as.numeric( determinant(b$d2K, logarithm=TRUE)$modulus )
  if( !log ) spa <- exp(spa)
  
  outList <- list("llk" = drop(spa), "mix" = mix)
  
  if(deriv)
  {
    grad <- .ecgf(lambda = lambda, X = X, kum1 = preCov$mY, kum2 = preCov$COV, mix = mix, grad = 2, 
                  deriv = TRUE, addList = list("y" = y, "decay" = decay), mixMethod = mixMethod )
    outList$grad <- grad
  }
  
  return(outList)
  
})


#### Additional details
#
### [con_test] Convergence Test:
# All components of "dK" have to be within "tol" after rescaling using the variances. 
# Tests with badly scaled systems suggests that the scaling is helpful. 
# Tests also say that scaling using the whole covariance matrix, is _NOT_ a good idea.