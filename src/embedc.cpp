#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix embedc(NumericVector x, int lags) {
  unsigned long int n = x.size();
  unsigned long int i = 0;
  unsigned long int j = 0;
  unsigned long int n_lags = n-lags;
  unsigned int lag_plus = lags + 1;
  NumericMatrix res(n_lags, lag_plus);

  for (j = 0; j < lag_plus; j++)
  {
    for (i = 0; i < n_lags; i++)
    {
      res(i, j) = x(lags + i - j);
    }
  }

  return (res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

 // /*** R
 // (y <- rnorm(10))
 // embedc(y, 1)
 // embedc(y, 2)
 // */
