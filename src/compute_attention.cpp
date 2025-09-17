#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
using namespace Rcpp;

// Helper function to compute median
double compute_median(NumericVector x) {
  int n = x.size();
  NumericVector x_sorted = clone(x);
  std::sort(x_sorted.begin(), x_sorted.end());

  if (n % 2 == 0) {
    return (x_sorted[n/2 - 1] + x_sorted[n/2]) / 2.0;
  } else {
    return x_sorted[n/2];
  }
}

// Helper function to compute standard deviation
double compute_sd(NumericVector x) {
  double mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
  double sq_sum = 0.0;
  for (int i = 0; i < x.size(); i++) {
    sq_sum += (x[i] - mean) * (x[i] - mean);
  }
  return std::sqrt(sq_sum / (x.size() - 1));
}

// 1. Fixed Cosine Attention (proper window-based similarity)
// [[Rcpp::export]]
NumericMatrix cosine_attention_cpp(NumericVector series, int window_size = 3) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    for (int j = 0; j <= t; j++) {
      double similarity = 0.0;

      // Define windows around t and j
      int start_t = std::max(0, t - window_size + 1);
      int start_j = std::max(0, j - window_size + 1);
      int end_t = t;
      int end_j = j;

      int length = std::min(end_t - start_t + 1, end_j - start_j + 1);

      if (length > 0) {
        double dot_product = 0.0;
        double norm1 = 0.0;
        double norm2 = 0.0;

        for (int k = 0; k < length; k++) {
          double val1 = series[start_t + k];
          double val2 = series[start_j + k];
          dot_product += val1 * val2;
          norm1 += val1 * val1;
          norm2 += val2 * val2;
        }

        if (norm1 > 1e-10 && norm2 > 1e-10) {
          similarity = dot_product / (std::sqrt(norm1) * std::sqrt(norm2));
        }
      }

      attention_matrix(t, j) = similarity;
      sum_weights += similarity;
    }

    // Normalize
    if (sum_weights > 1e-10) {
      for (int j = 0; j <= t; j++) {
        attention_matrix(t, j) /= sum_weights;
      }
    } else {
      // Fallback to uniform weights
      for (int j = 0; j <= t; j++) {
        attention_matrix(t, j) = 1.0 / (t + 1);
      }
    }
  }

  return attention_matrix;
}

// 2. Exponential Decay Attention (most robust for time series)
// [[Rcpp::export]]
NumericMatrix exponential_attention_cpp(NumericVector series, double decay_factor = 5.0) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    for (int j = 0; j <= t; j++) {
      double distance = t - j;
      double weight = std::exp(-distance / decay_factor);
      attention_matrix(t, j) = weight;
      sum_weights += weight;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_weights;
    }
  }

  return attention_matrix;
}

// 3. Fixed Dot Product Attention (using normalized values)
// [[Rcpp::export]]
NumericMatrix dot_product_attention_cpp(NumericVector series) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  // Normalize series to prevent scale issues
  double series_mean = std::accumulate(series.begin(), series.end(), 0.0) / n;
  double series_sd = compute_sd(series);
  NumericVector normalized_series(n);
  for (int i = 0; i < n; i++) {
    normalized_series[i] = (series_sd > 1e-10) ? (series[i] - series_mean) / series_sd : 0.0;
  }

  for (int t = 0; t < n; t++) {
    double max_score = -1e10;
    double sum_exp = 0.0;

    // Compute similarity scores and find max for numerical stability
    for (int j = 0; j <= t; j++) {
      double score = normalized_series[j] * normalized_series[t];
      if (score > max_score) max_score = score;
      attention_matrix(t, j) = score;
    }

    // Apply softmax
    for (int j = 0; j <= t; j++) {
      double exp_score = std::exp(attention_matrix(t, j) - max_score);
      attention_matrix(t, j) = exp_score;
      sum_exp += exp_score;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_exp;
    }
  }

  return attention_matrix;
}

// 4. Scaled Dot Product Attention (fixed with proper normalization)
// [[Rcpp::export]]
NumericMatrix scaled_dot_product_attention_cpp(NumericVector series, double temperature = 1.0) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  // Normalize series
  double series_mean = std::accumulate(series.begin(), series.end(), 0.0) / n;
  double series_sd = compute_sd(series);

  for (int t = 0; t < n; t++) {
    double max_score = -1e10;
    double sum_exp = 0.0;

    // Compute scaled similarity scores
    for (int j = 0; j <= t; j++) {
      double norm_j = (series_sd > 1e-10) ? (series[j] - series_mean) / series_sd : 0.0;
      double norm_t = (series_sd > 1e-10) ? (series[t] - series_mean) / series_sd : 0.0;
      double score = (norm_j * norm_t) / temperature;
      if (score > max_score) max_score = score;
      attention_matrix(t, j) = score;
    }

    // Apply softmax
    for (int j = 0; j <= t; j++) {
      double exp_score = std::exp(attention_matrix(t, j) - max_score);
      attention_matrix(t, j) = exp_score;
      sum_exp += exp_score;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_exp;
    }
  }

  return attention_matrix;
}

// 5. Gaussian Attention
// [[Rcpp::export]]
NumericMatrix gaussian_attention_cpp(NumericVector series, double sigma = 1.0) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    for (int j = 0; j <= t; j++) {
      double distance = std::abs(t - j);
      double weight = std::exp(-(distance * distance) / (2.0 * sigma * sigma));
      attention_matrix(t, j) = weight;
      sum_weights += weight;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_weights;
    }
  }

  return attention_matrix;
}

// 6. Linear Attention
// [[Rcpp::export]]
NumericMatrix linear_attention_cpp(NumericVector series) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    // Linear weights favoring recent observations
    for (int j = 0; j <= t; j++) {
      double weight = j + 1.0;  // Linear increasing weight
      attention_matrix(t, j) = weight;
      sum_weights += weight;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_weights;
    }
  }

  return attention_matrix;
}

// 7. Value-based Attention
// [[Rcpp::export]]
NumericMatrix value_based_attention_cpp(NumericVector series, double sensitivity = 1.0) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  // Compute robust scale (MAD - median absolute deviation)
  double median_val = compute_median(series);
  NumericVector abs_deviations(n);
  for (int i = 0; i < n; i++) {
    abs_deviations[i] = std::abs(series[i] - median_val);
  }
  double mad_scale = compute_median(abs_deviations);
  if (mad_scale < 1e-10) mad_scale = 1.0;

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    for (int j = 0; j <= t; j++) {
      double value_diff = std::abs(series[j] - series[t]) / mad_scale;
      double weight = std::exp(-value_diff * sensitivity);
      attention_matrix(t, j) = weight;
      sum_weights += weight;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_weights;
    }
  }

  return attention_matrix;
}

// 8. Hybrid Attention
// [[Rcpp::export]]
NumericMatrix hybrid_attention_cpp(NumericVector series, double time_decay = 5.0, double value_sensitivity = 1.0) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  // Robust scaling for values
  double median_val = compute_median(series);
  NumericVector abs_deviations(n);
  for (int i = 0; i < n; i++) {
    abs_deviations[i] = std::abs(series[i] - median_val);
  }
  double mad_scale = compute_median(abs_deviations);
  if (mad_scale < 1e-10) mad_scale = 1.0;

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    for (int j = 0; j <= t; j++) {
      double time_weight = std::exp(-(t - j) / time_decay);
      double value_diff = std::abs(series[j] - series[t]) / mad_scale;
      double value_weight = std::exp(-value_diff * value_sensitivity);
      double weight = time_weight * value_weight;
      attention_matrix(t, j) = weight;
      sum_weights += weight;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_weights;
    }
  }

  return attention_matrix;
}

// 9. Parametric Attention
// [[Rcpp::export]]
NumericMatrix parametric_attention_cpp(NumericVector series, double alpha = 0.5, double beta = 0.5) {
  int n = series.length();
  NumericMatrix attention_matrix(n, n);

  // Robust scaling
  double median_val = compute_median(series);
  NumericVector abs_deviations(n);
  for (int i = 0; i < n; i++) {
    abs_deviations[i] = std::abs(series[i] - median_val);
  }
  double mad_scale = compute_median(abs_deviations);
  if (mad_scale < 1e-10) mad_scale = 1.0;

  for (int t = 0; t < n; t++) {
    double sum_weights = 0.0;

    for (int j = 0; j <= t; j++) {
      double time_component = std::exp(-alpha * (t - j));
      double value_component = std::exp(-beta * std::abs(series[j] - series[t]) / mad_scale);
      double weight = time_component * value_component;
      attention_matrix(t, j) = weight;
      sum_weights += weight;
    }

    // Normalize
    for (int j = 0; j <= t; j++) {
      attention_matrix(t, j) /= sum_weights;
    }
  }

  return attention_matrix;
}

// Universal context vector computation
// [[Rcpp::export]]
NumericVector compute_context_vectors_cpp(NumericVector series, NumericMatrix attention_weights) {
  int n = series.length();
  NumericVector context_vectors(n);

  for (int t = 0; t < n; t++) {
    double weighted_sum = 0.0;

    for (int j = 0; j <= t; j++) {
      weighted_sum += attention_weights(t, j) * series[j];
    }

    context_vectors[t] = weighted_sum;
  }

  return context_vectors;
}

// Helper function to compare all mechanisms
// [[Rcpp::export]]
List compare_attention_mechanisms(NumericVector series,
                                  int window_size = 3,
                                  double decay_factor = 5.0,
                                  double temperature = 1.0,
                                  double sigma = 1.0,
                                  double sensitivity = 1.0) {
  List results;

  results["cosine"] = cosine_attention_cpp(series, window_size);
  results["exponential"] = exponential_attention_cpp(series, decay_factor);
  results["dot_product"] = dot_product_attention_cpp(series);
  results["scaled_dot_product"] = scaled_dot_product_attention_cpp(series, temperature);
  results["gaussian"] = gaussian_attention_cpp(series, sigma);
  results["linear"] = linear_attention_cpp(series);
  results["value_based"] = value_based_attention_cpp(series, sensitivity);
  results["hybrid"] = hybrid_attention_cpp(series, decay_factor, sensitivity);
  results["parametric"] = parametric_attention_cpp(series);

  return results;
}
