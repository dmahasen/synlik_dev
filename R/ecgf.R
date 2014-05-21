#####
#' Empirical cumulant generating function
#' @description Calculates the empirical cumulant generating function (CGF) and its derivatives
#'               given a sample of n d-dimentional vectors
#'
#' @param lambda Point at which the CGF is evaluated (d-dimensional vector).
#' @param X (n by d) matrix containing the data.
#' @param grad If grad == 0 only the value of the CGF is returned, 
#'             if grad == 1 also its first derivative wrt lambda 
#'             and if grad == 2 also the second derivarive wrt lambda.
#' @param mix Mixture of empirical and normal CGF to use (if 1 only empirical CGF is used).
#'            Default value is 0.9.
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{K} }{The value of the empirical CGF at lambda;}
#'         \item{ \code{dK} }{The value of the gradient empirical CGF wrt lambda at lambda;}
#'         \item{ \code{d2K} }{The value of the hessian of the empirical CGF wrt lambda at lambda;}
#'         }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#' @examples 
#' X <- matrix(rnorm(2 * 1e3), 1e3, 2)
#' ecgf(lambda = c(0, 0), X = X) 
#' @export
#'
ecgf <- function(lambda, X, grad = 0, mix = 0.9) {
  ## X[i,j] is ith rep of jth variable. Evaluate observed KGF 
  ## and its derivs w.r.t. lambda, without overflow...
  
  .ecgf(lambda = lambda, 
        X = X, 
        kum1 = colMeans(X), 
        kum2 = .robCov(t(X), alpha2 = 4, beta2 = 1.25)$COV, 
        grad = grad, 
        mix = mix )
  
}


#####
#' Empirical cumulant generating function
#' @description Calculates the empirical cumulant generating function (CGF) and its derivatives
#'               given a sample of n d-dimentional vectors
#'
#' @param lambda Point at which the CGF is evaluated (d-dimensional vector).
#' @param X (n by d) matrix containing the data.
#' @param kum1 Mean vector of the data.
#' @param kum2 Covariance matrix of the data.
#' @param grad If grad == 0 only the value of the CGF is returned, 
#'             if grad == 1 also its first derivative wrt lambda 
#'             and if grad == 2 also the second derivarive wrt lambda.
#' @param deriv If TRUE the gradient of the empitical CGF wrt y (and at y) is returned.
#'              Otherwise the values of the empirical CGF (and possibly of its derivatives wrt
#'              lambda) at lambda is returned.
#' @param mix Mixture of empirical and normal CGF to use (if 1 only empirical CGF is used).
#' @param addList = list of additional (optional) arguments: 
#'         \itemize{
#'         \item{ \code{invCOV} }{The inverse of kum2;}
#'         \item{ \code{y} }{The point at which the underlying empirical saddlepoint is evaluated;}
#'         \item{ \code{grad} }{The decay rate of the saddlepoint. See ?dsaddle for details;}
#'         }
#' @return If deriv == FALSE a list with entries:
#'         \itemize{
#'         \item{ \code{K} }{The value of the empirical CGF at lambda;}
#'         \item{ \code{dK} }{The value of the gradient empirical CGF wrt lambda at lambda;}
#'         \item{ \code{d2K} }{The value of the hessian of the empirical CGF wrt lambda at lambda;}
#'         }
#'         If deriv == TRUE the gradient of the empitical CGF wrt y (and at y) is returned. 
#'         This is used to calculate the gradient of the underlying empirical saddlepoint density
#'         at y. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#'


.ecgf <- cmpfun( function(lambda, X, kum1, kum2, mix, grad, deriv = FALSE, mixMethod = "mse", addList = list(NaN)) {
  ## X[i,j] is ith rep of jth variable. Evaluate observed KGF 
  ## and its derivs w.r.t. lambda, without overflow...
  
  if(grad > 2 && is.nan(addList[[1]])) stop("if you want the gradient you have to specify the additional arguments")
  
  if(!is.vector(lambda)) lambda <- as.vector(lambda) 
  if (!is.matrix(X)) X <- matrix(X, length(X), 1)
  n <- nrow(X)
  d <- ncol(X)
  
  ret <- .Call("ecgfCpp",
               lambda_ = lambda, 
               X_ = X, 
               mix_ = mix, 
               grad_ = grad, 
               kum1_ = kum1, 
               kum2_ = kum2,
               PACKAGE = "synlik")
  
  if(deriv)
  {
    lx <- lambda * t(X) 
    elx <- drop( exp( colSums(lx) - max(lx)) ) # exp(lambda'x_i - alpha) vector
    
    tmp_K <- ret$tmp_K
    tmp_dK <- drop( ret$tmp_dK )
    tmp_d2K <- ret$tmp_d2K
    
    switch(mixMethod,
           "gaus" = { dGam <- drop( (mix * addList$decay / sqrt(d)) * solve(kum2, kum1 - addList$y) )  },
           
           "mse" = { normRes <- solve(kum2, addList$y - kum1)
                     scaledRSS <- crossprod(addList$y - kum1, normRes) / sqrt(d)
                     sadVar <- exp(scaledRSS)
                     normVar <- scaledRSS + 0.5 * scaledRSS^2 + 1
                     d_sadVar <- 2 * sadVar * normRes / sqrt(d)
                     d_normVar <- 2 * normRes * ( 1 + scaledRSS ) / sqrt(d) 
                     dGam <- drop( d_normVar * sadVar - normVar * d_sadVar ) / (sadVar^2)
                     
                     if(normVar / sadVar < 1){
                       dGam <- addList$decay * ( (normVar / sadVar) ^ (addList$decay-1) ) * dGam 
                     } else {
                       dGam <- dGam * 0.0
                     }
           },
           stop("\"mixMethod\" must be either \"gaus\" or \"mse\"") 
    )
    
    if( addList$decay == Inf ) dGam <- numeric(d) 
    
    D <- diag( diag(ret$d2K) ^ -0.5, nrow = d, ncol = d)
    d2KQR <- qr( D %*% ret$d2K %*% D, tol = 0 )
    
    dLambda <- D %*% qr.solve(d2KQR, D %*% (diag(1, d) - tcrossprod(tmp_dK - kum1 - kum2%*%lambda, dGam)), tol = 0)
    
    d3K <- array(NA, c(d, d, d) )
    
    for(ff in 1:d)
    {
      A <- crossprod(X, (elx*X)*X[ , ff] ) / sum(elx)
      B <- -( tmp_d2K + tcrossprod(tmp_dK, tmp_dK) ) * tmp_dK[ff] 
      C <- - tcrossprod(tmp_d2K[ , ff], tmp_dK) - t( tcrossprod(tmp_d2K[ , ff], tmp_dK) )
      d3K[ , , ff] <- A + B + C 
    }  
    
    d2Kdy <- matrix(0, d, d)
    spaGrad <- numeric(d)
    for( ff in 1:d )
    {
      d2Kdy <- d2Kdy * 0 
      for(zz in 1:d)
      {
        d2Kdy <- d2Kdy + d3K[ , , zz] * dLambda[zz, ff] # Should it be [ff, zz] ??
      }
      
      spaGrad[ff] <- - mix * 0.5 * .Trace( D %*% qr.solve(d2KQR, D %*% d2Kdy, tol = 0) ) 
    }
    
    spaGrad <- spaGrad - 
      lambda + 
      dLambda %*% (ret$dK - addList$y) -  
      0.5 * .Trace( D %*% qr.solve(d2KQR, D %*% ( tmp_d2K - kum2 ), tol = 0) ) * dGam + 
      dGam * drop( tmp_K - crossprod(kum1, lambda) - 0.5 * crossprod(lambda, kum2%*%lambda) )
    
    return(spaGrad)
    
  } else return( ret[c("K", "dK", "d2K")] )
  
})