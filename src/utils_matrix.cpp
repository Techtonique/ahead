#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix create_lags_cpp(NumericVector x, int k)
{
  unsigned long int n = x.size();
  unsigned long int n_k = n-k;

  int ncol = k + 1;
  NumericMatrix res(n_k, ncol);

    for(int i = 0; i < ncol; i++) {
      for(int p = 0; p < n_k; p++)
      {
        res(p, i) = x(p + i);
      }
    }
    return res;
}

// [[Rcpp::export]]
List create_train_inputs_cpp(NumericMatrix x, int k)
{
  unsigned long int n = x.nrow();
  unsigned long int p = x.ncol();
  unsigned long int n_k = n-k;

  unsigned long int j, l;

  NumericMatrix y(n_k, p);
  NumericMatrix regressors(n_k, k*p);
  NumericMatrix Y_x(n_k, k);

  for(j = 0; j < p; j++) {
    Y_x = create_lags_cpp(x(_, j), k);
    y(_, j) = Y_x(_, 0);
      for(l = 1; l < k; l++) {
        regressors(_, j*k + (l-1)) = Y_x(_, l);
      }
      regressors(_, j*k + (k-1)) = Y_x(_, l);
  }

  return Rcpp::List::create(Rcpp::Named("y") = y,
                            Rcpp::Named("regressors") = regressors);

}

// [[Rcpp::export]]
NumericMatrix reformat_cpp(NumericMatrix x, unsigned long int n_k)
{
  unsigned long int n = x.nrow();
  unsigned long int p = x.ncol();
  NumericMatrix res(1, n_k*p);

  if (n_k >= n) {
    ::Rf_error("you must have n_k < x.nrow()");
  }

  for(int j = 0; j < p; j++) {
    for(int i = 0; i < n_k; i++) {
     res(0, j*n_k + i) = x(i, j);
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix rev_matrix_cpp(NumericMatrix x)
{
  unsigned long int n = x.nrow();
  unsigned long int p = x.ncol();
  NumericMatrix res(n, p);

int j;
 for(int k = 0; k < p; k++) {
   j = 0;
   for(int i = (n-1); i >= 0; i--) {
     res(j , k) = x(i , k);
     j++;
   }
 }
 return res;
}

// [[Rcpp::export]]
NumericMatrix rbind_vecmat_cpp(NumericVector y, NumericMatrix x)
{
  unsigned long int x_nrow = x.nrow();
  unsigned long int x_ncol = x.ncol();
  unsigned long int y_ncol = y.size();
  if (x_ncol != y_ncol) {
    ::Rf_error("you must have x.ncol() == y.size()");
  }
  unsigned long int res_nrow = x_nrow + 1;
  NumericMatrix res(res_nrow, x_ncol);

  for(int i = 0; i < x_nrow; i++) {
    res(i + 1 , _) = x(i , _);
  }

  for(int j = 0; j < x_ncol; j++) {
  res(0, j) = y(j);
  }

  return res;
}
