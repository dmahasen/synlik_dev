############################
##### Function to estimate gradient and hessian of synthetic likelihood
############################
# This is Simon's local regression to get gradient and hessian with one modification:
# the second regression (the one involving the residuals) is linear not quadratic. 
# This is suggested by the article "Local Polynomial Variance-Function Estimation" 
# since the variance shouldn't need non linear terms locally. 
# 
#' Estimate gradient and hessian of the synthetic likelihood.
#' 
#' @param param (numeric) vector containing the current values of the model's parameters.
#' @param nsim   (integer) number of simulated statistics.
#' @param covariance (matrix) covariance matrix used to simulate the parameters.
#' @param addRegr (logical). 
#'           If FALSE the statistics calculated by object@@summaries will be used (SL approach);
#'           if TRUE the simulated parameters will be regressed on the statistics and the 
#'           fitted values of the paramaters given the _observed_ statistics will be used as statistics
#'           (referred to as SL+ approach).
#' @param constr (named list) used to impose constraints on the parameters.
#'           Composed of 3 elements:
#'           [["indexes"]] = (numeric integers) indexes of the elements to check;
#'           [["upper"]]   = (numeric) upper bounds for the elements in "indexes";
#'           [["lower"]]   = (numeric) lower bounds for the elements in "indexes".
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to object@@simulator and object@@summaries.
#'            In general I would avoid using it and including in those two function everything they need.
#'
#' @return  a list containing: ["gradient"] = (numeric) estimated gradient of the synthetic log-likelihood at currPar;
#'                          ["hessian"]  = (matrix) estimated hessian of the synthetic log-likelihood at currPar;  
#'                          ["llk"]      = (scalar) estimated value of the synthetic log-likelihood at currPar.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                         
#' @export
#' 
synGrad_P <- cmpfun(function(param, nsim, covariance, 
                           simulator,d,obsData,summaryStat,uRN,
                           addRegr = TRUE, fixVar = TRUE, nBoot = 100, 
                           constr = list(), 
                           multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                           tolVar = 10 * .Machine$double.eps, ...)
{
  
  theData <- obsData
  
  nPar <- length(param)
  
  param <- unname(param)
  covariance <- unname(covariance)
  
 
  
  simulParams <- .paramsSimulator(theMean = param, covar = covariance, nsim = nsim, constr = constr)
  
  ##MAHSEN
  if(is.null(simulParams) | !is.numeric(simulParams))
  {
    stop(paste('simulParams is not numeric - theta = ', param, " obsY = ", obsData, sep=""))
  }
  
  
  simulData  <-  simulator(simulParams,nsim,d,uRN) 
  
  simulStats <- summaryStat(simulData,d)
    #.simulate.synlik(object, param = simulParams, nsim = nsim, stats = TRUE, clean = FALSE, verbose = FALSE, ...)  
  
  
  
  # Clean the simulated statistics from NaNs
  clean <- .clean(X = simulStats, verbose = TRUE)
  if(clean$nBanned > 0){
    simulStats <- clean$cleanX
    simulParams <- simulParams[-clean$banned, ]
  }
  rm(clean)
  
  #### Begin  MAHASEN 1/11/2018
  
  infStat <- is.infinite(rowSums(simulStats))
  
  simulStats <- simulStats[!infStat,]
  simulParams <- simulParams[!infStat,]
  
  numInfStat <- sum(infStat)
  if(numInfStat > 0.75 * nsim)
  {
    stop( paste(numInfStat, "out of", nsim, "statistics vectors", "contain INF and will not be used") )
  }
  
  

  if(!is.matrix(simulParams))
  {
    simulParams <- matrix(simulParams,ncol = 1)
  }
  
  ### End  MAHASEN 1/11/2018
  
  # Remove outliers that lay more than 6 MADs away from the median
  meds <- colMedians( simulStats )
  mads <- colMads( simulStats )
  outl <- colAnys( (abs((t(simulStats) - meds)) - 6 * mads) > 0 )
  if( any(outl) ){
    simulStats <- simulStats[ !outl, ]
    simulParams <- simulParams[ !outl, ]
  }
  
  # MAHASEN - to avoid errors in fitting regression model of summary statistics with parameters 
  #check summary statitics with variance zero. 
  badSStats <- which(diag(cov(simulStats))==0)
  
  if(length(badSStats)>0)
  {
    simulStats <- simulStats[,-badSStats]
  }
  
  nGood <- nrow(simulStats)
  
  # Adding the latest component to the mixture
  stored <- list("X" = simulParams, 
                 "y" = simulStats, 
                 "mu" = param, 
                 "sigma" = covariance)
  
  # if addRegr == TRUE we use SL+ and we regress parameters on statistics 
  if(addRegr){ 
    simulStats <- cbind(1, simulStats)
    nStats <- nPar
    
    qrx <- qr(simulStats, tol = 0)
    # matrix of linear regs coeffs for all summary statistics
    fearn_beta <- t( qr.coef(qrx, simulParams) )
    fearn_beta[is.na(fearn_beta)] <- 0 # Putting to zero the regression coefficients that are NAs
    
    simulStats <- tcrossprod(fearn_beta, simulStats)
    
   # summaries <- object@summaries
   # obserStats <- if( !is.null(summaries) ) drop( summaries(x = theData, extraArgs = object@extraArgs, ...) ) else drop( theData )
    obserStats <- summaryStat(theData,d)
    
    # MAHASEN
    if(length(badSStats)>0)
    {
      obserStats <- obserStats[-badSStats] 
    }
    obserStats <- fearn_beta %*% c(1, obserStats)
  } else{
    simulStats <- t(simulStats)
    nStats <- nrow(simulStats)
    
    #summaries <- object@summaries
   # obserStats <- if( !is.null(summaries) ) drop( summaries(x = theData, extraArgs = object@extraArgs, ...) ) else drop( theData )
    obserStats <- summaryStat(theData,d)
  }
  
  out <- .bootGrad(simulParams = simulParams, simulStats = simulStats, obserStats = obserStats, 
                   param = param, nBoot = nBoot, tolVar = tolVar, fixVar = fixVar)
  
  return( c(out, "stored" = list(stored)) )
  
})


# Boostrapps the simulated statistics and parameters to obtain the most
# positive definite Hessian
.bootGrad <- function(simulParams, simulStats, obserStats, param, nBoot, tolVar, fixVar)
{
  
  nPar <- length(param)
  nGood <- ncol(simulStats)
  
  # Create matrix X1 for quadratic regression of summary statistic (simulStats) 
  # on the parameters of each replicate (simulParams) 
  # Matrix X2 will be use to regress linearly the residuals of the first regression on the parameters
  
  X1 <- .quadModMat(simulParams, center = param, deriv = TRUE) # Model matrix for first regression 
  
  X2 <- X1[ , 1:(nPar+1)]                                      # Model matrix for second regression
  
  # Getting gradient and hessian using the full data
  out <- .synGrad(X1 = X1, X2 = X2, simulStats = simulStats, obserStats = obserStats, tolVar = tolVar, fixVar = fixVar)
  hessian <- out[["hessian"]]
  eig <- eigen(out[["hessian"]], only.values = TRUE)$values
  cond <- ifelse(all(eig < 0), max(abs(eig)) / min(abs(eig)), Inf)
  
  # Bootstrap to get a more negative definite hessian
  if(nBoot)
  {
    for(iBo in 1:nBoot)
    {
      # Sampling the model matrix
      index <- sample(1:nGood, nGood, replace = TRUE)
      bootX1 <- X1[index, ]
      bootX2 <- X2[index, ]
      bootStats <- simulStats[ , index, drop = FALSE]
      
      tmp <- .synGrad(X1 = bootX1, X2 = bootX2, simulStats = bootStats, 
                      obserStats = obserStats, tolVar = tolVar, fixVar = fixVar)
      tmpHess <- tmp[["hessian"]]
      
      eig <- eigen(tmpHess, only.values = TRUE)$values
      
      # Checking if the new hessian is better conditioned
      tmpCond <- ifelse(all(eig < 0), max(abs(eig)) / min(abs(eig)), Inf)
      
      if(tmpCond < cond){
        out[["hessian"]] <- tmpHess
        cond <- tmpCond
      }
      
    }
  }
  
  # If hessian is still not negative definit we make the diagonal entries all negative
  # if( !is.finite(cond) ) out[["hessian"]] <- out[["hessian"]] - diag( max(diag(out[["hessian"]])) + 1e-6, nPar)
  
  return( out )
  
}


# Given simulated parameters and statistics calculates gradient and hessian of the objective.
.synGrad <- function(X1, X2, simulStats, obserStats, tolVar, fixVar)
{
  
  nPar <- ncol(X2) - 1
  nStats <- nrow(simulStats)
  nSim <- ncol(simulStats)
  
  # Calculate the coefficients of the 1st regression
  qrx1 <- qr(X1, tol = 0)
  beta <- t( qr.coef(qrx1, t(simulStats)) )
  beta[ is.na(beta) ] <- 0
  
  ########################################################
  ######### Regressing entries of the covariance matrix on parameters 
  ########################################################
  
  # WARNING: it is better to simulated very close to the current value of parameters,
  # otherwise the quadratic model is biased, we get mean(res) != 0 and SIGMA_HAT is not invertible!
  
  # Obtaining residuals (res) in the two ways previously described #res is a (nStats X nSIm) matrix 
  res <- simulStats - tcrossprod(beta, X1) 
  
  # Down-weighting and discarding the most extreme residuals
  fixCov <- robCov(res)
  
  # Matrix of expected SIGMA, first derivative of sigma wrt the parameters
  # and second derivatives. (Example: SIGMA12 = (D^2 SIGMA)/(D theta1 D theta2))     
  SIGMA_HAT <- fixCov$COV
  firstSIG <- array(0, c(nStats, nStats, nPar))
  
  # Now I regress each vector D_vett[ , , i] on the parameters (using the same LINEAR regression)
  # The intercepts are the expected elements of the covariance matrix. I move along the upper triangular
  # part of the matrix. 
  # N.B. Since this is a LINEAR regression this code is not the same as Simon's notes. 
  # Calculate the derivatives of the variances (the diagonal), not the whole covariance matrix.
  
  if( !fixVar )
  {
    
    resSd <- fixCov$sd
    
    varCoef <- expRegr(y = t((res / resSd)^2), X = X2[ , -1, drop = F], weights = fixCov$weights) 
    
    varCoef[1, ] <- varCoef[1, ] + 2 * log(resSd)
    
    diag(SIGMA_HAT) <- exp( varCoef[1, ] )
    
    COR <- cov2cor(fixCov$COV)
    
    for(ii in 1:nStats)
      for(jj in ii:nStats)
      {
        firstSIG[ii, jj, ] <- 
          firstSIG[jj, ii, ]  <- COR[ii, jj] * exp(0.5*(varCoef[1, ii]+varCoef[1, jj])) * 0.5*(varCoef[-1, ii] + varCoef[-1, jj])
      }
    
    tmp <- sqrt( diag(SIGMA_HAT) )
    
    SIGMA_HAT <- diag(tmp, nStats) %*% COR %*% diag(tmp, nStats)
  }
  
  # Identify statistics with very low variance (diagonal elements of SIGMA_HAT)
  # and removing them and their coefficients
  lowVar <- which( diag(SIGMA_HAT) < tolVar )
  if( length(lowVar) ) {
    if( addRegr ) stop( paste("The variance of one of the parameters is  <", tolVar, ""))
    nStats <- nStats - length(lowVar)
    if(nStats == 0) stop( paste("All the chosen statistics have variance <", tolVar) )
    beta <- beta[-lowVar, ]
    obserStats <- obserStats[-lowVar]
    firstSIG <- firstSIG[-lowVar, -lowVar, ] 
    SIGMA_HAT <- SIGMA_HAT[-lowVar, -lowVar] 
    warning( paste("There are", length(lowVar), "statistics with variance <", tolVar, "they were removed.") )
  }
  
  # Get the (scaled) QR decomposition of SIGMA_HAT, that later I'll use to 
  # solve several linear systems (so I avoid inverting it).
  D <- diag(diag(SIGMA_HAT)^-0.5, nStats)
  
  sigQR <- qr( D%*%SIGMA_HAT%*%D, tol = 0 )
  
  #############################################################################
  ###### Calculate gradient and Hessian of the log-likelhood wrt the parameters
  #############################################################################
  
  ######### Calculating the gradient of the synthetic likelihood
  # DmuDth1 = the first derivatives of mu wrt the parameters 
  # DmuDth2 = the second derivatives of mu wrt the parameters
  # currStat = the expected value of the statistic at the current position
  # obsRes = statistics(observed_path) - curr_stat
  # sigByRes = SIGMA_HAT^-1 %*% (y - currStat)
  DmuDth1 <- beta[ , 2:(nPar+1), drop = FALSE]
  DmuDth2 <- beta[ , (nPar+2):ncol(X1), drop = FALSE]
  currStat <- beta[ , 1]
  obsRes <- drop( obserStats - currStat )
  sigByRes <- drop(D %*% qr.solve(sigQR, D %*% obsRes, tol = 0))
  
  firstDeriv <- numeric(nPar)
  for(kk in 1:nPar)
  {
    firstDeriv[kk] <- 0.5 * (2 * crossprod(DmuDth1[ , kk], sigByRes)    +     
                               
                               crossprod(sigByRes, firstSIG[ , , kk]%*%sigByRes) -
                               
                               .Trace(D %*% qr.solve(sigQR, D %*% firstSIG[ , , kk], tol = 0)) )
  }
  
  
  ######### Calculating the Hessian of the synthetic likelihood
  # Different from Simon's document because the second regression is linear, and hence
  # all the second derivatives of the covariance matrix wrt the paramets have been put to zero.
  # I fill the hessian moving on the upper triangle and I use the switch to select the correct
  # matrix of the first derivatives of the covariance wrt the parameters. 
  
  secondDeriv <- matrix(0, nPar, nPar)
  
  if(nPar == 1){
    
    indexes <- matrix(1, 1, 1)
    
  }else{
    
    # Create matrix of indexes to manage the second derivarives (stored in DmuDth2)
    indexes <- diag(seq(1:nPar))
    entries <- seq(nPar + 1, nPar + factorial(nPar)/(factorial(2)*factorial(nPar-2)))
    zz <- 1
    for(jj in 1:nPar){
      indexes[jj, -(1:jj)] <- entries[zz:(zz + nPar - jj - 1)]
      zz <- zz + nPar - jj
    }
  }
  
  
  for(kk in 1:nPar) 
    for(ff in kk:nPar)
    {
      zz <- indexes[kk, ff]
      
      DsigDthK <- firstSIG[ , , kk]
      DsigDthL <- firstSIG[ , , ff]
      
      secondDeriv[kk, ff] <- secondDeriv[ff, kk] <-
        0.5 * (  
          2 * crossprod(DmuDth2[ , zz], sigByRes) - 
            2 * crossprod(DmuDth1[ , kk], D %*% qr.solve(sigQR, D%*%(DsigDthL%*%sigByRes), tol = 0)) -
            2 * crossprod(DmuDth1[ , kk], D %*% qr.solve(sigQR, D%*%DmuDth1[ , ff], tol = 0)) - 
            2 * crossprod(DmuDth1[ , ff], D %*% qr.solve(sigQR, D%*%(DsigDthK%*%sigByRes), tol = 0)) -  
            2 * crossprod(sigByRes, DsigDthL%*%(D %*% qr.solve(sigQR, D%*%(DsigDthK%*%sigByRes), tol = 0)) ) +
            .Trace(D %*% qr.solve(sigQR, D%*%DsigDthL, tol = 0) %*% (D %*% qr.solve(sigQR, D%*%DsigDthK, tol = 0)) )
        )  
    }
  
  llk <- - 0.5 * crossprod(obsRes, sigByRes) - 0.5 * log( abs( prod( diag(qr.R(sigQR)) / diag(D^2) ) ) )
  
  return( list("gradient" = firstDeriv, "hessian" = secondDeriv, "llk" = llk) )
  
}