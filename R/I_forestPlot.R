
# Simple forest plot:
# mu is the mean vector
# sigma the vector of standard error
# xlabs are the labels on the x axis

.forestPlot <- function(mu, sigma, xlabels, ...)
{
  n <- length(mu)
  
  sigma[ is.na(sigma) ] <- 0.0
  
  upbar <- mu + sigma
  lowbar <- mu - sigma
  
  ylims <- c( min(lowbar), max(upbar) ) 
  
  plot(1:n, mu, xaxt = "n", pch = 19, ylim = ylims, ...)
  axis(1, at=1:n, labels = xlabels)
  
  for(ii in 1:n)
    segments(x0 = ii, 
             y0 = lowbar[ii],
             x1 = ii,
             y1 = upbar[ii], 
             lwd = 2)
  
  abline(h = min(mu), lty = 2, lwd = 2) 
  
  return( invisible(NULL) )
}


# synlik:::.forestPlot(1:3, 3:1, xlabels = c("a", "b", "c"))
