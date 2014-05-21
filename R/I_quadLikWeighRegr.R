
.quadLikWeighRegr <- function(llk, parMat, nreps, boot, quant = 0, ...)
{ 
  # If "parMat" is a 3d array, we need to transform it into a matrix.
  if( is.matrix(llk) ) llk <- as.vector(llk)
  if( length(dim(parMat)) == 3 ){
    
    parMat <- alply(parMat, 3, function(input) if(is.matrix(input)) input else as.matrix(input))
    parMat <- do.call("rbind", parMat)

  }
  
  npar <- ncol(parMat)
  linPar <- 1 : (npar+1)
  
  # Non finite loglikelihoods are excluded
  good <- is.finite(llk)
  parMat <- parMat[good, , drop = FALSE]
  llk <- llk[ good ]
  
  # Discarding values corresponding to positions that are too far from the center.
  if( quant )
  {    
    mu <- apply(parMat, 2, median)
    Sig <- cov(parMat)
    
    dist <- maha(X = parMat, mu = mu, sigma = Sig)
    
    good <- ( dist < qchisq(quant, npar) )
    parMat <- parMat[good, , drop = FALSE]
    llk <- llk[ good ]
  } 
  
  nval <- length(llk)
  
  # Creating model matrix and weights
  X <- .quadModMat(parMat)
  
  # Checking the sample size
  if( nval < npar * 10 ) warning(paste("Sample size = ", nval, "but you have", npar, "parameters. \n",
                                        "Maybe better getting more likelihood evaluations."))
  
  # Regressing estimated log-likelihoods of parameters
  fit <- rlm(llk ~ -1 + X, maxit = 100)
  
  # Extracting quadratic and linear coefficients
  bestCoef <- coef(fit)
  lcoef <- bestCoef[ linPar ]
  qcoef <- bestCoef[ -linPar ]
  ncoef <- length(bestCoef)
  nqcoef <- length(qcoef)
  
  # Extract hessian and its eigen-values
  hess <- - .extractHessian(qcoef, npar)
  eig <- eigen(hess, only.values = TRUE)$values
  
  ###### Start of PD correction
  # If hessian is not PD, bootstrap or simulate using asymptotic covariance to get it PD (hopefully).
  # Our final estimate is the positive definite hessian with the lowest conditioning number.
  if( any(eig < 0) )
  {
    message("The estimated Hessian of the log-lik is not negative definite, I will try to bootstrap it.")
    
    # Resampled coefficients will be stored here by row
    coefMat <- matrix(NA, nreps, ncoef)
    
    # Using boostrap
    if(boot)
    {
      for(ii in 1:nreps)
      {
        index <- sample(1:nval, nval, replace = TRUE)
        
        tmpX <- X[index, ]
        tmpllk  <- llk[index]
        
        coefMat[ii, ] <- lm.fit(x = tmpX, y = tmpllk)$coefficients
      }
    } else {
      # Using asymptotic covariance
      coefMat <- .rmvn(nreps, mu = allCoef, sigma = vcov(fit)) 
    }
    
    # Extract simulated hessian and their conditioning number
    cond <- aaply(coefMat[ , -linPar, drop = FALSE], 
                  1, 
                  function(inVett){
                    negH  <- - .extractHessian(inVett, npar)
                    eig <- eigen(negH, only.values = TRUE)$values
                    if( all(eig > 0) ){ 
                      return( max(eig) / min(eig) ) 
                    } else { 
                      return( NA ) 
                    } 
                  })
    
    if( all(is.na(cond)) ) stop("Cannot get a positive definite negative Hessian.")
    
    # Use hessian with lowest conditioning number
    bestCoef <- coefMat[which.min(cond), ]
    hess <- - .extractHessian(bestCoef[-linPar], npar)
  } 
  
  return( list("linear" = bestCoef[linPar], "quadratic" = - hess / 2, "hessian" = -hess) )
  
}

