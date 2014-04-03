
.extractHessian <- cmpfun(function(quadCoef, nPar)
{
  hessian <- matrix(NA, nPar, nPar)
  
  # Calculate indexes to store coefficients in the Hessian
  if(nPar == 1){
    
    indexes <- matrix(1, 1, 1)
    
  }else{
    
    # Create matrix of indexes to manage the second derivarives stored in beta
    indexes <- diag(seq(1:nPar))
    entries <- seq(nPar + 1, nPar + factorial(nPar)/(factorial(2)*factorial(nPar-2)))
    zz <- 1
    for(jj in 1:nPar){
      indexes[jj, -(1:jj)] <- entries[zz:(zz + nPar - jj - 1)]
      zz <- zz + nPar - jj
    }
  } 
  
  # Put derivatives into the hessian.
  # Diagonal multiplied by 2
  for(iRow in 1:nPar)
    for(iCol in iRow:nPar)
      hessian[iRow, iCol] <- hessian[iCol, iRow] <- quadCoef[ indexes[iRow, iCol] ] * ifelse(iRow == iCol, 2, 1) 
  
  return( hessian )
  
})