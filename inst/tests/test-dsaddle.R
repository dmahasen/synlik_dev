
context("dsaddle")

test_that( "Testing saddlepoint density mode-finding", {
  
  #########
  #### Univariate case
  #########
  
  # Test if it gets close to the mode
  x <- rgamma(10000, 2, 1)
  decay <-  0.2
  
  res_mse <- findMode(x, init = 8, decay = decay, sadTol = 1e-6, mixMethod = "mse")$mode
  res_gaus <- findMode(x, init = 8, decay = decay, sadTol = 1e-6, mixMethod = "gaus")$mode
  
  expect_less_than( abs(res_mse - 1.03), 1e-1 )
  expect_less_than( abs(res_gaus - 1.2), 1e-1 )
  
  # Testing gradient with finite differences
  delta <- 10^-6
  posit <- 6
  fd_mse <- - dsaddle(y = posit + delta, X = x, tol = 10^-6, decay = decay, mixMethod = "mse", log = T)$llk + 
    dsaddle(y = posit, X = x, tol = 10^-6, decay = decay, mixMethod = "mse", log = T)$llk
  fd_gaus <- - dsaddle(y = posit + delta, X = x, tol = 10^-6, decay = decay, mixMethod = "gaus", log = T)$llk + 
    dsaddle(y = posit, X = x, tol = 10^-6, decay = decay, mixMethod = "gaus", log = T)$llk
  
  sadGrad_mse <- -dsaddle(y = posit + delta/2, X = x, tol = 10^-6, decay = decay, mixMethod = "mse", deriv = TRUE)$grad
  sadGrad_gaus <- -dsaddle(y = posit + delta/2, X = x, tol = 10^-6, decay = decay, mixMethod = "gaus", deriv = TRUE)$grad
  
  expect_less_than( abs(fd_mse/delta - sadGrad_mse), 1e-3 )
  expect_less_than( abs(fd_gaus/delta - sadGrad_gaus), 1e-3 )
  
  
  #########
  #### Multivariate case
  #########
  library(mvnfast)
  # Test if it gets close to the mode
  decay <- 0.2
  dims <- 3
  A <- matrix(rnorm(dims^2), dims, dims)
  A <- t(A)%*%A
  myMu <- seq.int(0, dims-1)*10
  SIGMA <-  A
  X <- rmvn(n = 2000, mu = myMu, sigma = SIGMA)
  
  init <- colMeans(X) + rmvn(n = 1, mu = rep(0, dims), sigma = SIGMA)

  mode_mse <- findMode(X = X, init = drop(init), decay = decay, mixMethod = "mse")$mode
  mode_gaus <- findMode(X = X, init = drop(init), decay = decay, mixMethod = "gaus")$mode
  
  expect_true( all(abs(mode_mse - colMeans(X)) / sqrt(diag(SIGMA)) < 1e-1) )
  expect_true( all(abs(mode_gaus - colMeans(X)) / sqrt(diag(SIGMA)) < 1e-1) )
  
  # Testing gradient with finite differences
  dims <- 6
  decay <- 0.5
  A <- matrix(rnorm(dims^2), dims, dims)
  A <- t(A)%*%A
  myMu <- seq.int(1, dims)
  SIGMA <- A
  X <- rmvn(n = 10000, mu = myMu, sigma = SIGMA)
  delta <- 1e-6
  
  fd_mse <- fd_gaus <- sad_grad_mse <- sad_grad_gaus <- numeric(dims)
  for(index in 1:dims)
  {
    pert <- rep(0, dims); 
    pert[index] <- delta
    x = colMeans(X) + rmvn(n = 1, mu = rep(0, dims), sigma = SIGMA)
    
    aa <-  dsaddle(y = x + pert, X = X, tol = 10^-6, decay, log = TRUE, mixMethod = "mse")$llk - 
      dsaddle(y = x, X = X, tol = 10^-6, decay, log = TRUE, mixMethod = "mse")$llk
    fd_mse[index] <- aa/delta
    sad_grad_mse[index] <- dsaddle(y = x + pert/2, X = X, tol = 10^-6, decay = decay, deriv = TRUE, mixMethod = "mse")$grad[index]
    
    aa <-  dsaddle(y = x + pert, X = X, tol = 10^-6, decay, log = TRUE, mixMethod = "gaus")$llk - 
      dsaddle(y = x, X = X, tol = 10^-6, decay, log = TRUE, mixMethod = "gaus")$llk
    fd_gaus[index] <- aa/delta
    sad_grad_gaus[index] <- dsaddle(y = x + pert/2, X = X, tol = 10^-6, decay = decay, deriv = TRUE, mixMethod = "gaus")$grad[index]
  }
    
  expect_true( all( abs(sad_grad_mse - fd_mse) < 1e-3) )
  expect_true( all( abs(sad_grad_gaus - fd_gaus) < 1e-3) )
  
}
)


