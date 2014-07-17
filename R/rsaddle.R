########
# Simulating r.v. from a saddlepoint density
########


rsaddle <- function(n, X, decay, ncores = 1, ...)
{
  prop <- rmvn(n, colMeans(X), cov(X), ncores = ncores)
  
  w <- dsaddle(prop, X = X, decay = decay, ...)$llk / dmvn(prop, colMeans(X), cov(X), ncores = ncores)
    
  out <- prop[ sample(1:n, n, replace = TRUE, prob = w), ]
  
  return( out )
  
}