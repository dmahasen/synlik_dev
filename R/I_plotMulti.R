
# Plotting list of "smcmc" objects

.plotMulti <- function(inlist, chainSlot = "chains", truePar = NULL, ...)
{
  nchains <- length(inlist)
  npar <- ncol( slot(inlist[[1]], chainSlot) )
  
  parNames <- names(inlist[[1]]@param)
  
  for(ipar in 1:npar)
  {
    tmp <- sapply(inlist, function(input) slot(input, chainSlot)[ , ipar])
    matplot(tmp, type = 'l', ylab = parNames[ipar], ...)
    if( !is.null(truePar) ) abline(h = truePar[ipar], col = 1, lwd = 3, ...)
  }
  
  tmpMean <- lapply(inlist, function(input) colMeans(slot(input, chainSlot)))
  tmpCov <- lapply(inlist, function(input) cov(slot(input, chainSlot)))
  
  mu <- Reduce("+", tmpMean) / nchains
  covar <- Reduce("+", tmpCov) / nchains
  
  return( invisible( list("mean" = mu, "covar" = covar) ) )
}