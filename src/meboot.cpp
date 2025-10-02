#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector meboot_approx(NumericVector p, int n, NumericVector z, 
                            NumericVector desintxb) {
  NumericVector q(n);
  int nz = z.size();
  
  for (int i = 0; i < n; i++) {
    double pi = p[i];
    bool found = false;
    // Check intervals for linear interpolation
    // Interval j: [(j)/n, (j+1)/n]
    for (int j = 0; j < nz - 1; j++) {
      double lower = (double)j / n;
      double upper = (double)(j + 1) / n;
      
      if (pi > lower && pi <= upper) {
        // Linear interpolation
        double prop = (pi - lower) / (upper - lower);
        q[i] = z[j] + prop * (z[j + 1] - z[j]);
        found = true;
        break;
      }
    }
    // If not found in main intervals, handle edge cases
    if (!found) {
      if (pi <= (1.0 / n)) {
        // Left tail - will be handled in main function
        q[i] = NA_REAL;
      } else if (pi >= ((double)(n - 1) / n)) {
        // Right tail - will be handled in main function
        q[i] = NA_REAL;
      } else {
        // Fallback: simple linear interpolation across all z
        int idx = (int)(pi * (nz - 1));
        idx = std::max(0, std::min(idx, nz - 2));
        double prop = pi * (nz - 1) - idx;
        q[i] = z[idx] + prop * (z[idx + 1] - z[idx]);
      }
    }
  }
  
  return q;
}

// [[Rcpp::export]]
NumericMatrix meboot_part_rcpp(NumericVector x, int reps, 
                               NumericVector z, double xmin, 
                               double xmax, NumericVector desintxb, 
                               bool reachbnd) {
  int n = x.size();
  int nz = z.size();
  NumericMatrix ensemble(n, reps);
  // Use R's RNG for reproducibility with set.seed()
  GetRNGstate();
  
  for (int j = 0; j < reps; j++) {
    NumericVector p(n);
    // Generate random numbers using R's RNG
    for (int i = 0; i < n; i++) {
      p[i] = R::runif(0.0, 1.0);
    }
    // Main interpolation
    NumericVector q = meboot_approx(p, n, z, desintxb);
    // Handle tails
    for (int i = 0; i < n; i++) {
      double pi = p[i];
      // Left tail
      if (pi <= (1.0 / n)) {
        double prop = pi * n;
        q[i] = xmin + prop * (z[0] - xmin);
        if (!reachbnd) {
          q[i] = q[i] + desintxb[0] - 0.5 * (z[0] + xmin);
        }
      }
      // Right tail  
      else if (pi >= ((double)(n - 1) / n)) {
        double prop = (pi - ((double)(n - 1) / n)) * n;
        q[i] = z[nz - 1] + prop * (xmax - z[nz - 1]);
        if (!reachbnd) {
          int last_idx = desintxb.size() - 1;
          q[i] = q[i] + desintxb[last_idx] - 0.5 * (z[nz - 1] + xmax);
        }
      }
    }
    
    ensemble(_, j) = q;
  }
  
  PutRNGstate();
  
  return ensemble;
}

// Helper function for expand.sd
// [[Rcpp::export]]
NumericMatrix expand_sd_rcpp(NumericVector x, NumericMatrix ensemble, 
                             double fiv = 1.0, double elbow = 0.95, 
                             bool force_clt = false) {
  int n = x.size();
  int reps = ensemble.ncol();
  NumericMatrix result(n, reps);
  // Calculate statistics for original series
  double xmed = median(x);
  // Materialize the expression into a vector before taking median
  NumericVector x_abs = abs(x - xmed);
  double xmad = median(x_abs);
  // Calculate mean and SD for original series (for force_clt option)
  double xmean = mean(x);
  double xsd = sd(x);
  for (int j = 0; j < reps; j++) {
    NumericVector y = ensemble(_, j);
    if (force_clt) {
      // Use mean and SD expansion (Central Limit Theorem)
      double ymean = mean(y);
      double ysd = sd(y);
      double ratio = (ysd > 1e-10) ? (xsd / ysd) : 1.0;
      for (int i = 0; i < n; i++) {
        result(i, j) = xmean + (y[i] - ymean) * ratio;
      }
    } else {
      // Use median and MAD expansion (more robust)
      double ymed = median(y);
      // Materialize the expression into a vector before taking median
      NumericVector y_abs = abs(y - ymed);
      double ymad = median(y_abs);
      // Calculate expansion factor using elbow parameter
      double ratio = (ymad > 1e-10) ? (xmad / ymad) : 1.0;
      double expand_factor = elbow + (1.0 - elbow) * ratio;
      for (int i = 0; i < n; i++) {
        result(i, j) = xmed + (y[i] - ymed) * expand_factor;
      }
    }
  }
  return result;
}

// Alternative version that computes medians more efficiently
// [[Rcpp::export]]
NumericMatrix expand_sd_rcpp_fast(NumericVector x, NumericMatrix ensemble, 
                                  double fiv = 1.0, double elbow = 0.95) {
  int n = x.size();
  int reps = ensemble.ncol();
  NumericMatrix result(n, reps);
  // Pre-calculate statistics for original series
  double xmed = median(x);
  // Materialize the expression before taking median
  NumericVector x_abs = abs(x - xmed);
  double xmad = median(x_abs);
  for (int j = 0; j < reps; j++) {
    NumericVector y = ensemble(_, j);
    double ymed = median(y);
    // Materialize the expression before taking median
    NumericVector y_abs = abs(y - ymed);
    double ymad = median(y_abs);
    // Calculate expansion factor
    double ratio = (ymad > 1e-10) ? (xmad / ymad) : 1.0;
    double expand_factor = elbow + (1.0 - elbow) * ratio;
    // Apply expansion
    for (int i = 0; i < n; i++) {
      result(i, j) = xmed + (y[i] - ymed) * expand_factor;
    }
  }
  return result;
}