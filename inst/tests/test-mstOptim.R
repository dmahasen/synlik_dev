context("mstOptim")

test_that("mstOptim works well with slik and selectDecay", {
  
  expSimul <- function(param, nsim, extraArgs, ...)
  {
    
    if( !all( c("nObs") %in% names(extraArgs) ) ) stop("extraArgs should contain nObs")
    nObs <- extraArgs$nObs
    
    if( !is.loaded("plyr") ) library("plyr")
    
    if( !is.matrix(param) ) param <- matrix(param, nsim, length(param), byrow = TRUE)
    
    d <- ncol(param)
    
    output <- alply(param, 
                    1,
                    function(inPar){
                      matrix(rexp(d*nObs, inPar), nObs, d, byrow = TRUE)
                    })
    
    return( output )
  }
  
  # Summaries
  expStats <- function(x, extraArgs, ...){
    ## obsData is a vector of observed path
    ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
    
    if(!is.list(x)) x <- list(x)
    
    X0 <- laply(x, colMeans)
    
    return(X0) 
  }

  # SL object  
  exp_sl <- new("synlik", 
                simulator = expSimul,
                summaries = expStats,
                param = c("alpha1" = 1, "alpha2" = 2), 
                extraArgs = list("nObs" = 200)
  )
  
  #### Simulate from the object
  param <- exp_sl@param
  exp_sl@data <- simulate(exp_sl, clean = FALSE, nsim = 1)[[1]]
  
  set.seed(4414)
  slow <- slik(exp_sl, 
               param  = param,
               nsim    = 200, saddle = TRUE,  
               controlSad = list("decay" = 1, "normalize" = TRUE, "method" = "IS", "nNorm" = 1e3, "fastInit" = FALSE) )
  
  set.seed(4414)
  fast <- slik(exp_sl, 
               param  = param,
               nsim    = 200, saddle = TRUE,
               controlSad = list("decay" = 1, "normalize" = TRUE, "method" = "IS", "nNorm" = 1e3, "fastInit" = TRUE) )
    
  expect_less_than( abs(slow - fast) / abs(slow), 1e-6 )
  
  set.seed(4414)
  cvSlow <- selectDecay(decay = c(0.05, 0.5, 2, 10), 
                        simulator = function(...) simulate(exp_sl, nsim = 1000, stats = TRUE), 
                        K = 2, 
                        normalize = T,
                        control = list("method" = "IS", "nNorm" = 100) )
  
  set.seed(4414)
  cvFast <- selectDecay(decay = c(0.05, 0.5, 2, 10), 
                        simulator = function(...) simulate(exp_sl, nsim = 1000, stats = TRUE), 
                        K = 2, 
                        normalize = T,
                        control = list("method" = "IS", "nNorm" = 100, "fastInit" = TRUE) )
    
  expect_less_than( max(abs(as.vector(cvSlow$negLogLik - cvFast$negLogLik)) / abs(as.vector(cvSlow$negLogLik))), 1e-6)
  
  
})