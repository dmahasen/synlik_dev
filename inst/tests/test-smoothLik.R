
context("smoothLik - coeff.sgmm")

test_that( "Additive regression (mgcv) is able to approximate the synthetic likelihood correctly", {
  
  ################################################################################################################
  ################################### FIRST TEST                               ###################################
  ################################################################################################################
  
  # E(S) = mu = beta * param
  # COV(S) = sigma = diag( exp( alpha + gamma * param )  )
  # N simulated parameters and statistics
  # d dimension of statistics and parameters
  
  library(mvnfast)
  library(plyr)
  
  N <- 2e4
  d <- 2
  
  paramTrue <- runif(d, -10, 10)
  simPar <- rmvn(N, paramTrue, diag(1:d, d))
  
  beta <- 1 : d
  alpha <- 1
  gamma <- 0.1 * rep(1, d)
  
  scaledPar <- scale(simPar)
  center <- attr(scaledPar, "scaled:center")
  scale <- attr(scaledPar, "scaled:scale")
  
  varFun <- function(inPar){ 
    if(d != 2) stop("This version works only with d == 2")
    
    stdDev <- sqrt( diag( drop(exp(alpha + diag(drop(scale(t(inPar), center, scale)), d) %*% gamma)), d) )
    
    return( stdDev %*% matrix(c(1, 0.9, 0.9, 1), 2, 2) %*% stdDev )
  }
  
#   varFun <- function(inPar){ 
#     diag( drop(exp(alpha + diag(drop(scale(t(inPar), center, scale)), d) %*% gamma)), d)
#   }
  
  covs <- alply(simPar, 1, varFun)
  mu <- t(t(simPar) * beta)
  simStat <- mu + laply(covs, function(inCov) rmvn(1, rep(0, d), inCov))
  
  statTrue <- beta * paramTrue
    
  normLik <- smoothLik(yobs = statTrue, y = simStat, param = simPar, multicore = T, ncores = 2)
  
  # Test estimation of gmm objective function
  dsl <- 50
  xSeq <- matrix(paramTrue, dsl, d, byrow = TRUE)
  xSeq <- xSeq*matrix(seq(0.5, 2, length = dsl), dsl, d)
  
  smo <- apply(xSeq, 1, normLik, fun = "gmm")
  truth <- apply(xSeq,
                 1,
                 function(input)
                   t(statTrue - beta * input) %*% solve(varFun(input)) %*% (statTrue - beta * input))
  
  plot(xSeq[ , 1], smo, type = 'l', main = "GMM objective (quadratic form)", ylab = "GMM")
  lines(xSeq[ , 1], truth, col = 2)
  legend("topleft", c("smoothed", "truth"), col = c(1, 2), lty = 1) 
  
  expect_less_than( mean(abs(smo - truth)), 1)
  
  # Test estimation of synthetic likelihood
  smo <- apply(xSeq, 1, normLik, fun = "sl")
  truth <- apply(xSeq,
                 1,
                 function(input) drop( dmvn(statTrue, beta * input, varFun(input), log = TRUE ) ))
  
  plot(xSeq[ , 1], smo, type = 'l', main = "Synthetic Likelihood", ylab = "SL")
  lines(xSeq[ , 1], truth, col = 2)
  legend("topright", c("smoothed", "truth"), col = c(1, 2), lty = 1) 
  
  expect_less_than( mean(abs(smo - truth)), 1)
  
  
  
  
  ################################################################################################################
  ################################### SECOND TEST                               ##################################
  ################################################################################################################
  
  #######
  # Synlik model with multivariate normal model, the statistics are the samples means
  ########
  
  # Simulator
  normSimul <- function(param, nsim, extraArgs, ...)
  {
    
    if( !all( c("nObs") %in% names(extraArgs) ) ) stop("extraArgs should contain nObs")
    nObs <- extraArgs$nObs
    covarFun <- extraArgs$covarFun
    
    if( !is.loaded("plyr") ) library("plyr")
    
    if( !is.matrix(param) ) param <- matrix(param, nsim, length(param), byrow = TRUE)
    
    d <- ncol(param)
    
    output <- alply(param, 
                    1,
                    function(inPar){
                      rmvn(nObs, inPar, covarFun(inPar))
                    })
    
    return( output )
  }
  
  # Statistics
  normStats <- function(x, extraArgs, ...){
    ## obsData is a vector of observed path
    ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
    
    if(!is.list(x)) x <- list(x)
    
    X0 <- laply(x, colMeans)
    
    return(X0) 
  }
  
  ######################
  # Starting test
  #######################
  
  d <- 3
  
  mu <- 1:d
  a <- runif(1, 0.01, 1)
  b <- runif(d, 0.1, 0.5)
  
  names(mu) <- sapply(1:d, function(input) paste("mu", input, sep = '') )
  
  covarFun <- function(inPar){ diag( exp(a + inPar * b) ) }
  
  norm_sl <- new("synlik", 
                 simulator = normSimul,
                 summaries = normStats,
                 param = mu, 
                 extraArgs = list("nObs" = 200, 
                                  "covarFun" =  covarFun)
  )
  
  #### Simulate from the object
  norm_sl@data <- simulate(norm_sl, clean = FALSE, nsim = 1)[[1]]
  
  nsim <- 2e4
  par <- synlik:::.paramsSimulator(norm_sl@param, 
                                   covarFun(norm_sl@param) * 2, 
                                   nsim, 
                                   constr = list("indexes" = 1:d, 
                                                 "upper" = rep(Inf, d), 
                                                 "lower" = rep(0, d)))
  
  stat <- norm_sl@simulator(par, nsim = nsim, extraArgs = norm_sl@extraArgs)
  stat <- norm_sl@summaries(x = stat, extraArgs = norm_sl@extraArgs)
  
  obsStat <- norm_sl@summaries(x = norm_sl@data, extraArgs = norm_sl@extraArgs)
  
  normLik <- smoothLik(yobs = obsStat, y = stat, param = par, multicore = TRUE, ncores = 2)
  
  dsl <- 50
  xSeq <- matrix(norm_sl@param, dsl, d, byrow = TRUE)
  xSeq <- xSeq + matrix(seq(-0.5, 0.5, length = dsl), dsl, d)
  
  smo <- apply(xSeq, 1, normLik, fun = "gmm")
  truth <- apply(xSeq,
                 1,
                 function(input)
                   maha(input, obsStat, covarFun(input) / norm_sl@extraArgs[["nObs"]]))
  plot(xSeq[ , 1], smo, type = 'l', main = "GMM objective (quadratic form)", ylab = "GMM")
  lines(xSeq[ , 1], truth, col = 2)
  legend("top", c("smoothed", "truth"), col = c(1, 2), lty = 1) 
  
  expect_less_than( mean(abs(smo - truth)), 1)
  
  smo <- apply(xSeq, 1, normLik, fun = "sl")
  truth <- apply(xSeq,
                 1,
                 function(input)
                   dmvn(input, obsStat, covarFun(input) / norm_sl@extraArgs[["nObs"]], d, log = TRUE))
  plot(xSeq[ , 1], smo, type = 'l', main = "GMM objective (quadratic form)", ylab = "GMM")
  lines(xSeq[ , 1], truth, col = 2)
  legend("topleft", c("smoothed", "truth"), col = c(1, 2), lty = 1) 
  
  expect_less_than( mean(abs(smo - truth)), 1)
  
})










