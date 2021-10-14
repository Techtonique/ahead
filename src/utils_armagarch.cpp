#include <Rcpp.h>
using namespace Rcpp;


// Forecasting innovations - core loop

// [[Rcpp::export]]
NumericVector forecast_innovs_loop_cpp(NumericVector eps,
                               NumericVector rts,
                               double eps_prev,
                               double omega, double alpha, double beta,
                               unsigned int df, unsigned int h)
{
  double sigma_prev, sigma_new;
  double eps_new;

  sigma_prev = omega/(1 - alpha - beta);
  sigma_new = 0;
  eps_new = 0;

  for(unsigned int i = 0; i < h; i++)
  {
    sigma_new = sqrt(omega + alpha*pow(eps_prev, 2) + beta*pow(sigma_prev, 2));
    eps_new = rts(i)*sigma_new/sqrt(df/(df-2));

    eps(i) = eps_new;

    eps_prev = eps_new;
    sigma_prev = sigma_new;
  }

  return(eps);
}


// [[Rcpp::export]]
NumericVector forecast_innovs_loop_cpp2(NumericVector eps,
                                       NumericVector rn,
                                       double eps_prev,
                                       double omega, double alpha, double beta,
                                       unsigned int df, unsigned int h)
{
  double sigma_prev, sigma_new;
  double eps_new;

  sigma_prev = omega/(1 - alpha - beta);
  sigma_new = 0;
  eps_new = 0;

  for(unsigned int i = 0; i < h; i++)
  {
    sigma_new = sqrt(omega + alpha*pow(eps_prev, 2) + beta*pow(sigma_prev, 2));
    eps_new = rn(i)*sigma_new;

    eps(i) = eps_new;

    eps_prev = eps_new;
    sigma_prev = sigma_new;
  }

  return(eps);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
