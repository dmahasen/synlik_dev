# Plotting list of "smcmc" objects

.plotChains <- function(ch, pnam, tp = NULL, ...)
{  
  nch <- length(ch)
  npar <- ncol( ch[[1]] )
  
  for(ipar in 1:npar)
  {
    tmp <- sapply(ch, function(input) input[ , ipar])
    matplot(tmp, type = 'l', ylab = pnam[ipar], ...)
    if( !is.null(tp) ) abline(h = tp[ipar], col = 1, lwd = 3, ...)
  }
  
  means <- lapply(ch, function(input) colMeans(input))
  covs <- lapply(ch, function(input) cov(input))
  
  mu <- Reduce("+", means) / nch
  covar <- Reduce("+", covs) / nch
  
  names(mu) <- pnam
  
  return( invisible( list("mean" = mu, "covar" = covar) ) )
}