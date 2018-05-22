#include "synlik.h"

SEXP crashModelCpp(SEXP days_, SEXP nSimul_, SEXP param_, SEXP randInit_, SEXP initVal_)
{
  using namespace Rcpp;
  
  int days = as<int>(days_);
  int nSimul = as<int>(nSimul_);
  NumericMatrix param = as<NumericMatrix>(param_);
  bool randInit = as<bool>(randInit_);
  double initVal = as<double>(initVal_);
  
  int nparam = param.ncol(); 
  bool multiparam = false;
  
  if(nparam != 4) stop("Wrong number of parameters");
  if(param.nrow() > 1) { multiparam = true; }
  if(multiparam == true && param.nrow() != nSimul) 
    stop("Number of parameters vectors is different from the number of simulations");
  
  NumericMatrix output( nSimul, days );
  RNGScope scope; // Declare RNGScope after the output in order to avoid a known Rcpp bug.
  
  double r = exp(param(0, 0));
  double K = exp(param(0, 1));
  double alpha = exp(param(0, 2));
  double beta = exp(param(0, 3));
  
  NumericVector initState(nSimul);
  if(randInit){ initState = rpois(nSimul, 10); } else { initState = initState + initVal; }
  
  NumericVector::iterator initIter = initState.begin();
  
  double currState;
  
  for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
  {
    if( multiparam == true )
    {
      r = exp(param(iRow, 0));
      K = exp(param(iRow, 1));
      alpha = exp(param(iRow, 2));
      beta = exp(param(iRow, 3));
    }
    
    currState = *initIter;
    output(iRow, 0) = currState;

    for(int iCol = 1; iCol < days; iCol++){
      currState>K ? currState=rbinom(1,currState,alpha)[0] : currState=rpois(1,currState*(1+r))[0];
      currState += rpois(1,beta)[0];
      output(iRow, iCol) = currState;
    }

  }
  
  return output;
  
}