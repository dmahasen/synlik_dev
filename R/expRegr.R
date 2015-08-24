 
# This fits the regression y ~ exp( alpha + X %*% beta ) * z^2, where z ~ N(0, sigma).
# First get initial estimates by OLS on regression model log(y) ~ alpha + X %*% beta + log(z^2)
# Then we do maximum likelihood on true model, where model log-likelihood = log( Chi-square density * Jacobian ).
# So the log likelihood would be:
#
# p(y) = ChiSq(z) / Jacob = ChiSq{y / exp(alpha + X %*% beta)} / exp(alpha + X %*% beta)
#
# But the Jacobian can be absorbed into the density, so that
#
# y ~ Gamma(shape = 1/2, scale = 2 * exp( alpha + X %*% beta ))

expRegr <- function(y, X, weights = NULL, tol = 6.0)
{
  if( !is.matrix(y) )
  { 
    y <- matrix(y, length(y), 1)
  }
  
  if( !is.data.frame(X) ) 
  {
    if( any(is.null(colnames(X))) )
    { 
      colnames(X) <- sapply(1:ncol(X), function(input) paste("x", input, sep = ''))
    }
    X <- as.data.frame(X)
  }
  
  d <- ncol(y)
  np <- ncol(X)
    
  logY <- log(y)

  # Create model formulas for lm() and glm()
  lmForm <- list()
  lmForm[[1]] <- paste("tmpY ~ 1", sep = '')
  lmForm <- append(lmForm, lapply(1:np, function(input) paste("+ x", input, sep = '')))
  glmForm <- lmForm
  glmForm[[1]] <- paste("y[ , ii] ~ 1", sep = '')
  
  # Transforming lists into formulas
  lmForm <- as.formula( do.call("paste", lmForm) )
  glmForm <- as.formula( do.call("paste", glmForm) )
    
  fitCoef <- matrix(NA, np+1, d)
  for(ii in 1:d) # Loop on columns of y
  {
    
    ### FIRST FIT with OLS on log scale to get good initial values
    # Discard outliers
    tmpY <- logY[ , ii]
    good <- abs(tmpY - median(tmpY)) / mad(tmpY)  < tol
    tmpY <- tmpY[ good ]
    tmpX <- X[good, , drop = F]
    
    tmpCoef <- lm(lmForm, data = tmpX)$coef
    
    # Assuming the errors are chisq, by working on log scale we get a bias in the intercept equal to -1.270322...
    tmpCoef[1] <- tmpCoef[1] + 1.270322
    
    ### SECOND FIT with gamma glm with log link   
    fitCoef[ , ii] <- gam(glmForm, data = X, family = "Gamma"(link='log'), start = tmpCoef, weights = weights)$coef

  }
  
  if( !is.matrix(fitCoef) ){ fitCoef <- matrix(fitCoef, length(fitCoef), 1) }
        
  return(fitCoef)
  
}