
.plotIter <- function(A, trans = NULL, type = "line", addPoints = NULL, addPlot = NULL, ...)
{
  if( !(type %in% c("line", "hist")) ) stop("type should be either line or hist")
  
  if( is.vector(A) ) A <- array(A, dim = c(length(a), 1, 1), dimnames = list(NULL, names(A), NULL))
  if( is.matrix(A) ) A <- array(A, dim = c(dim(A), 1), dimnames = list(NULL, colnames(A), NULL))
  
  parNames <- dimnames(A)[[2]]

  if( dim(A)[1] > 0 ) 
  {
    nPar <- dim(A)[2]
    panelDim <- min( ceiling(sqrt(nPar)), 3 )
    par(mfrow = c(panelDim, panelDim))
    
    # Transform each matrix if needed
    A <- alply(A, 3, function(input) .transMatrix(input, trans))
    
    # Don't invert the following lines
    Amat <- do.call("rbind", A)
    A <- do.call("abind", c(A, "along" = 3))
    
    if( !is.null(addPoints) ){
      xpoints <- addPoints$x
      ypoints <- addPoints$y
      ypoints <- .transMatrix(ypoints, trans)
    }
    
    counter <- 1
    # Plot each column
    for(nam in parNames){
      
      # Either lines
      if(type == "line")
      {
        if( is.null(addPoints) )
        {
        
        ii <- which(parNames == nam)
        matplot(1:dim(A)[1], A[ , ii, ], type = 'l', main = nam,
             ylab = nam, xlab = "Iteration", ...)
        } else {
          plot(xpoints, ypoints[ , nam], main = nam,
               ylab = nam, xlab = "Iteration",
               ylim =  c(min(c(Amat[ , nam], ypoints[ , nam])), max(c(Amat[ , nam], ypoints[, nam]))), ...)
          lines(1:nrow(Amat), Amat[ , nam], col = 2, lwd = 2)
        }
        
      } else {
        # Or histograms
        hist(Amat[ , nam],  main = nam, ylab = nam, xlab = nam, ...)
      }
      
      if( !is.null(addPlot) ) get(addPlot)(nam, ...)
      
      if( !(counter %% (panelDim^2)) && (counter != nPar) ) readline(prompt = "Press <Enter> to see the next plot...") 
      counter <- counter + 1
    }
  }
  
  return(NULL)   
}