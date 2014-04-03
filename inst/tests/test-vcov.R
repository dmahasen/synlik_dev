
context("vcov")

test_that( "Quadratic regression is able to recover the covariance, maximal log-likelihood and MLE", {
  
  N <- 1000
  mu <- - 1
  beta <- 1:3
  
  Quad <- - matrix(c(3, 1, 0, 1, 3, 0, 0, 0, 1), 3, 3, byrow = T)
  Hess <- 2 * Quad
  Cov <- solve(-Hess)
  
  objFun <- function(x)
  {
    x <- drop(x)
    mu + t(beta) %*% x + t(x) %*% Quad %*% x
  }
  
  truePar <- - 0.5 * solve(Quad, beta) 
  
  y <- matrix(NA, N, 2)
  x <- array(rnorm(2 * N * 3), c(N, 3, 2))
  y[ , 1] <- apply(x[ , , 1], 1, objFun) + rnorm(N, 0, 1e-6) 
  y[ , 2] <- apply(x[ , , 2], 1, objFun) + rnorm(N, 0, 1e-6) 

  mcmcObj <- list("initPar" = matrix(c("a" = 1, "b" = 2, "c" = 3), 2, 3, byrow = T),
                  "propCov" = diag( c(0.2, 0.1, 0.2) ), 
                  "llkStore" = y,
                  "parStore" = x)
  
  # Covariance test
  estCov <- vcov.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than(sum(abs(estCov - Cov)), 1e-3)
  
  # Loglik test
  llkMax <- logLik.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than( abs(llkMax - objFun(truePar)) , 1e-3)
  
  # Estimates test
  estPar <- coef.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than( sum( abs(estPar - truePar) ) , 1e-3)
  
  ##############
  # Tests with one parameter blocked
  ##############
  
  x[ , 2, ] <- truePar[2]
  Cov[c(1, 3), c(1, 3)] <- solve(-Hess[c(1, 3), c(1, 3)])
  Cov[2, ] <- Cov[ , 2] <- 0
  y[ , 1] <- apply(x[ , , 1], 1, objFun) + rnorm(N, 0, 1e-6) 
  y[ , 2] <- apply(x[ , , 2], 1, objFun) + rnorm(N, 0, 1e-6) 
  
  mcmcObj <- list("initPar" = matrix(c("a" = 1, "b" = truePar[2], "c" = 3), 2, 3, byrow = T),
                  "propCov" = diag( c(0.2, 0, 0.2) ), 
                  "llkStore" = y,
                  "parStore" = x)
  
  # Covariance test
  estCov <- vcov.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than(sum(abs(estCov - Cov)), 1e-3)
  
  # Loglik test
  llkMax <- logLik.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than( abs(llkMax - objFun(truePar)) , 1e-3)
  
  # Estimates test
  estPar <- coef.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than( sum( abs(estPar - truePar) ) , 1e-3)
  
  ##############
  # Tests with two parameters blocked
  ##############
  
  x[ , 2, ] <- truePar[2]
  x[ , 3, ] <- truePar[3]
  Cov <- Cov * 0
  Cov[1, 1] <- solve(-Hess[1, 1])
  y[ , 1] <- apply(x[ , , 1], 1, objFun) + rnorm(N, 0, 1e-6) 
  y[ , 2] <- apply(x[ , , 2], 1, objFun) + rnorm(N, 0, 1e-6) 
  
  mcmcObj <- list("initPar" = matrix(c("a" = 1, "b" = truePar[2], "c" = truePar[3]), 2, 3, byrow = T),
                  "propCov" = diag( c(0.2, 0, 0) ), 
                  "llkStore" = y,
                  "parStore" = x)
  
  # Covariance test
  estCov <- vcov.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than(sum(abs(estCov - Cov)), 1e-3)
  
  # Loglik test
  llkMax <- logLik.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than( abs(llkMax - objFun(truePar)) , 1e-3)
  
  # Estimates test
  estPar <- coef.mcmc(object = mcmcObj, method = "mle", nreps = 1000, boot = TRUE)
  
  expect_less_than( sum( abs(estPar - truePar) ) , 1e-3)

})


