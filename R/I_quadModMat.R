
# Creates model matrix for quadratic multivariate regression 

.quadModMat <- function(X, center = FALSE, scale = FALSE, deriv = FALSE)
{
  #nPar <- ncol(X) # when number of paramters = 1 , X is a vector 
  
  X <- scale(X, center = center, scale = scale)
  nPar <- ncol(X)
  
  M <- cbind(1, X, X ^ 2 / ifelse(deriv, 2, 1))
     
  # Creating all mixed terms
  if(nPar > 1){
    
    comb <- combn(nPar, 2)
    
    tmp <- lapply(1:ncol(comb), function(jj) X[ , comb[1, jj]] * X[ , comb[2, jj]])
    
    M <- cbind(M, do.call("cbind", tmp))
    
  }
  
  return( M ) 
} 