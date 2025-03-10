#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_attention_cpp(NumericVector series) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);
  
  // Pre-compute series norm once since we'll use the full series
  double series_norm = 0.0;
  for (int i = 0; i < n; i++) {
    series_norm += series[i] * series[i];
  }
  series_norm = sqrt(series_norm);
  
  // For each time point t, compute attention weights for all points j <= t
  for (int t = 0; t < n; t++) {
    for (int j = 0; j <= t; j++) {
      double dot_product = 0.0;
      double past_norm = 0.0;
      
      // Only consider values up to point j
      for (int k = 0; k <= j; k++) {
        dot_product += series[k] * series[k];
        past_norm += series[k] * series[k];
      }
      past_norm = sqrt(past_norm);
      
      // Calculate similarity
      double similarity = 0.0;
      if (past_norm > 0.0) {
        similarity = dot_product / (series_norm * past_norm);
      }
      
      attention_matrix(t, j) = similarity;
    }
    
    // Normalize weights for current time point
    double sum_weights = 0.0;
    for (int j = 0; j <= t; j++) {
      sum_weights += attention_matrix(t, j);
    }
    
    if (sum_weights > 0.0) {
      for (int j = 0; j <= t; j++) {
        attention_matrix(t, j) /= sum_weights;
      }
    }
  }
  
  return attention_matrix;
}

// [[Rcpp::export]]
NumericVector compute_context_vectors_cpp(NumericVector series, 
                                        NumericMatrix attention_weights) {
  int n = series.length();
  NumericVector context_vectors(n);
  
  // For each time point t, compute weighted sum of past values
  for (int t = 0; t < n; t++) {
    double weighted_sum = 0.0;
    // Only use information up to time t (causal attention)
    for (int j = 0; j <= t; j++) {
      weighted_sum += attention_weights(t, j) * series[j];
    }
    context_vectors[t] = weighted_sum;
  }
  
  return context_vectors;
} 