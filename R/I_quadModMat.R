
# Creates model matrix for quadratic multivariate regression 

.quadModMat <- function(X)
{
  nPar <- ncol(X)
  
  M <- cbind(1, X, X ^ 2)
     
  # Creating all mixed terms
  if(nPar > 1){
    
    comb <- combn(nPar, 2)
    
    tmp <- lapply(1:ncol(comb), function(jj) X[ , comb[1, jj]] * X[ , comb[2, jj]])
    
    M <- cbind(M, do.call("cbind", tmp))
    
  }
  
  return( M ) 
} 