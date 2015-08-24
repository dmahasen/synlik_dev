context("expRegr")

test_that( "Testing non-linear exponential regression", {
  
  library(mvnfast)
  N <- 2e4
  d <- 3
  ny <- 2
  
  X <- rmvn(N, 1:d, diag(d:1, d))
  y <- c()
  
  b <- matrix(runif(ny*d, 0.1, 0.5), d, ny)
  a <- runif(ny, -1, 1)
  
  for(ii in 1:ny){
    y <- cbind(y, matrix(exp(a[ii] + X %*% b[ , ii]) * rnorm(N)^2) )
  }
  
  colnames(X) <- sapply(1:d, function(input) paste("x", input, sep = ''))
  X <- as.data.frame(X)
  
  fitCoef <- expRegr(y = y, X = X, tol = 5)
  
  xSeq <- matrix(0, 50, nrow(b), byrow = TRUE)
  xSeq <- xSeq + matrix( seq(-2, 2, length = 50), 50, nrow(b))  
  colnames(xSeq) <- names(X)
  
  pred <- 2 * exp( cbind(1, xSeq) %*% fitCoef )
  
  for(ii in 1:ny)
  {   
    truth <- 2 * exp(a[ii] + xSeq %*% b[ , ii])
    plot(pred[ , ii], type = 'l', main = "Exponential Regression", xlab = "x", ylab = "Predicted")
    lines(truth, col = 2)
    legend("topleft", c("smoothed", "truth"), col = c(1, 2), lty = 1) 
    expect_less_than( mean(abs(pred[ , ii] - truth)), 0.1 * mean(truth))
  }
  
#   mean(y[ , 1])
#   gam(y[ , 1] / mean(y[ , 1]) ~ 1 + x1 + x2 + x3, data = X, family = "Gamma"(link='log'))$coef
#   
#   gam(y[ , 1] ~ 1 + x1 + x2 + x3, data = X, family = "Gamma"(link='log') )$coef
  
})


# x <- rchisq(1e5, df = 1)
# 
# sigma = 100
# y <- sigma * x
# xSeq <- seq(0, max(y), length.out = 10000)
# hist(y, prob = T, breaks = 100)
# lines(xSeq, dgamma(xSeq, shape = 0.5, scale = 2 * sigma))
# 
# lines(xSeq, dchisq(xSeq / sigma, df = 1) / sigma, col = 2)
# 
# lines(xSeq, dchisq(xSeq, 1), col = 2)

