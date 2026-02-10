#' Generate synthetic time series via model-based residual bootstrap
#'
#' @param y Time series data (ts or mts object)
#' @param model_type One of "auto", "auto.arima", "tslm", "ets", "arima", "custom"
#' @param model_func Custom model function (required if model_type = "custom")
#' @param formula Formula for tslm (default: trend + season)
#' @param order ARIMA order (p,d,q) for model_type = "arima"
#' @param seasonal List with order and period for seasonal ARIMA
#' @param length_out Length of synthetic series to generate (default = length of y)
#' @param n_sim Number of synthetic series to generate
#' @param seed Random seed for reproducibility (optional)
#' @param ... Extra arguments passed to the modelling function
#'
#' @return A list with synthetic series and model information
#' @export
generate_synthetic_ts <- function(y,
                                  model_type = "auto",
                                  model_func = NULL,
                                  formula = NULL,
                                  order = NULL,
                                  seasonal = list(order = c(0,0,0), period = NA),
                                  length_out = NULL,
                                  n_sim = 1000,
                                  seed = NULL,
                                  ...) {
  
  # Input validation
  if (!inherits(y, "ts")) {
    stop("y must be a ts or mts object")
  }
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Store original time series properties
  ts_start <- start(y)
  ts_frequency <- frequency(y)
  is_multivariate <- is.mts(y)
  original_length <- if (is_multivariate) nrow(y) else length(y)
  
  # Set output length to match input if not specified
  if (is.null(length_out)) {
    length_out <- original_length
  }
  
  # Choose model function based on model_type and data type
  final_model_func <- switch(
    model_type,
    "auto" = {
      if (is_multivariate) {
        function(x, ...) {
          if (is.null(formula)) {
            tslm(x ~ trend + season, ...)
          } else {
            tslm(formula, data = x, ...)
          }
        }
      } else {
        forecast::auto.arima
      }
    },
    "auto.arima" = {
      if (is_multivariate) {
        stop("auto.arima only supports univariate time series. Use model_type = 'tslm' for multivariate.")
      }
      forecast::auto.arima
    },
    "tslm" = {
      function(x, ...) {
        if (is.null(formula)) {
          tslm(x ~ trend + season, ...)
        } else {
          tslm(formula, data = x, ...)
        }
      }
    },
    "ets" = {
      if (is_multivariate) {
        stop("ets only supports univariate time series.")
      }
      forecast::ets
    },
    "arima" = {
      if (is_multivariate) {
        stop("Arima only supports univariate time series. Use model_type = 'tslm' for multivariate.")
      }
      if (is.null(order)) {
        stop("order must be provided when model_type = 'arima'")
      }
      function(x, ...) Arima(x, order = order, seasonal = seasonal, ...)
    },
    "custom" = {
      if (is.null(model_func)) {
        stop("model_func must be provided when model_type = 'custom'")
      }
      model_func
    },
    stop("Invalid model_type. Choose from: 'auto', 'auto.arima', 'tslm', 'ets', 'arima', 'custom'")
  )
  
  # Fit model with na.exclude to maintain time alignment
  fit <- final_model_func(y, ...)
  
  # Extract fitted values and residuals
  fitted_vals <- fitted(fit)
  residuals_obj <- residuals(fit)
  
  # Generate synthetic series
  synthetic_series_list <- vector("list", n_sim)
  
  if (is_multivariate) {
    # MULTIVARIATE CASE - handle each series separately
    n_series <- ncol(y)
    series_names <- colnames(y)
    if (is.null(series_names)) {
      series_names <- paste0("Series", 1:n_series)
    }
    
    # Pre-compute residuals for each series
    residual_list <- vector("list", n_series)
    fitted_list <- vector("list", n_series)
    
    for (j in 1:n_series) {
      resid_j <- residuals_obj[, j]
      fitted_j <- fitted_vals[, j]
      valid_idx <- !is.na(resid_j) & !is.na(fitted_j)
      
      if (sum(valid_idx) == 0) {
        stop("No valid residuals for series ", series_names[j])
      }
      residual_list[[j]] <- resid_j[valid_idx]
      fitted_list[[j]] <- fitted_j[valid_idx]
    }
    
    for (sim in 1:n_sim) {
      sim_series <- matrix(NA, nrow = length_out, ncol = n_series)
      
      for (j in 1:n_series) {
        # Sample from fitted values + bootstrapped residuals
        fitted_clean <- fitted_list[[j]]
        resid_clean <- residual_list[[j]]
        
        # If we need longer series than original, sample with replacement
        if (length_out > length(fitted_clean)) {
          boot_indices <- sample(seq_along(fitted_clean), size = length_out, replace = TRUE)
          base_values <- fitted_clean[boot_indices]
          boot_residuals <- resid_clean[boot_indices]
        } else {
          boot_indices <- sample(seq_along(fitted_clean), size = length_out, replace = FALSE)
          base_values <- fitted_clean[boot_indices]
          boot_residuals <- resid_clean[boot_indices]
        }
        
        sim_series[, j] <- base_values + boot_residuals
      }
      
      # Convert to proper mts object with same time properties
      sim_ts <- ts(sim_series, start = ts_start, frequency = ts_frequency)
      colnames(sim_ts) <- series_names
      synthetic_series_list[[sim]] <- sim_ts
    }
    
  } else {
    # UNIVARIATE CASE
    # Clean fitted values and residuals
    valid_idx <- !is.na(fitted_vals) & !is.na(residuals_obj)
    fitted_clean <- fitted_vals[valid_idx]
    resid_clean <- residuals_obj[valid_idx]
    
    if (length(fitted_clean) == 0) {
      stop("No valid fitted values or residuals from the model")
    }
    
    for (sim in 1:n_sim) {
      # Sample from fitted values + bootstrapped residuals
      if (length_out > length(fitted_clean)) {
        # For longer series, sample with replacement
        boot_indices <- sample(seq_along(fitted_clean), size = length_out, replace = TRUE)
        base_values <- fitted_clean[boot_indices]
        boot_residuals <- resid_clean[boot_indices]
      } else {
        # For same or shorter length, sample without replacement
        boot_indices <- sample(seq_along(fitted_clean), size = length_out, replace = FALSE)
        base_values <- fitted_clean[boot_indices]
        boot_residuals <- resid_clean[boot_indices]
      }
      
      synthetic_ts <- base_values + boot_residuals
      synthetic_ts <- ts(synthetic_ts, start = ts_start, frequency = ts_frequency)
      synthetic_series_list[[sim]] <- synthetic_ts
    }
  }
  
  # Return results
  result <- list(
    synthetic_series = synthetic_series_list,
    original_series = y,
    fitted_values = fitted_vals,
    residuals = residuals_obj,
    model = fit,
    n_sim = n_sim,
    length_out = length_out,
    is_multivariate = is_multivariate,
    ts_properties = list(
      start = ts_start,
      frequency = ts_frequency
    ),
    model_info = list(
      model_type = model_type,
      formula = if (!is.null(formula)) formula else NULL,
      order = if (!is.null(order)) order else NULL
    )
  )
  
  class(result) <- c("synthetic_ts", "list")
  return(result)
}

#' Plot method for synthetic_ts objects
#'
#' @param x synthetic_ts object
#' @param type Type of plot: "series", "residuals", "acf", "density", "qq"
#' @param which_sim Which simulation to plot (for type = "series")
#' @param series_index Which series to plot (for multivariate data)
#' @param ... Additional arguments passed to plot functions
#'
#' @export
plot.synthetic_ts <- function(x, type = "series", which_sim = 1, series_index = 1, ...) {
  if (!inherits(x, "synthetic_ts")) {
    stop("Object must be of class 'synthetic_ts'")
  }
  
  # Set up plotting parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  switch(type,
         "series" = {
           .plot_series_comparison(x, which_sim, series_index, ...)
         },
         "residuals" = {
           .plot_residuals(x, series_index, ...)
         },
         "acf" = {
           .plot_acf_comparison(x, which_sim, series_index, ...)
         },
         "density" = {
           .plot_density_comparison(x, which_sim, series_index, ...)
         },
         "qq" = {
           .plot_qq_comparison(x, which_sim, series_index, ...)
         },
         stop("Invalid type. Choose from: 'series', 'residuals', 'acf', 'density', 'qq'")
  )
}

#' Summary method for synthetic_ts objects
#'
#' @param object synthetic_ts object
#' @param ... Additional arguments
#'
#' @export
summary.synthetic_ts <- function(object, ...) {
  if (!inherits(object, "synthetic_ts")) {
    stop("Object must be of class 'synthetic_ts'")
  }
  
  # Perform hypothesis tests
  tests <- .perform_adequacy_tests(object)
  
  # Create summary structure
  summary_obj <- list(
    basic_info = .get_basic_info(object),
    model_info = .get_model_info(object),
    statistical_tests = tests,
    goodness_of_fit = .get_goodness_of_fit(object)
  )
  
  class(summary_obj) <- "summary.synthetic_ts"
  return(summary_obj)
}

#' Print method for summary.synthetic_ts objects
#'
#' @param x summary.synthetic_ts object
#' @param ... Additional arguments
#'
#' @export
print.summary.synthetic_ts <- function(x, ...) {
  cat("=== SYNTHETIC TIME SERIES SUMMARY ===\n\n")
  
  # Basic information
  cat("BASIC INFORMATION:\n")
  cat("  Original series length:", x$basic_info$original_length, "\n")
  cat("  Synthetic series length:", x$basic_info$synthetic_length, "\n")
  cat("  Number of simulations:", x$basic_info$n_sim, "\n")
  cat("  Frequency:", x$basic_info$frequency, "\n")
  cat("  Multivariate:", x$basic_info$is_multivariate, "\n")
  if (x$basic_info$is_multivariate) {
    cat("  Number of series:", x$basic_info$n_series, "\n")
  }
  cat("\n")
  
  # Model information
  cat("MODEL INFORMATION:\n")
  cat("  Model type:", x$model_info$model_type, "\n")
  if (!is.null(x$model_info$formula)) {
    cat("  Formula:", deparse(x$model_info$formula), "\n")
  }
  if (!is.null(x$model_info$order)) {
    cat("  ARIMA order:", paste(x$model_info$order, collapse = ","), "\n")
  }
  cat("\n")
  
  # Goodness of fit
  cat("GOODNESS OF FIT (Original vs Synthetic):\n")
  cat(sprintf("  Mean Absolute Error: %.4f\n", x$goodness_of_fit$mae))
  cat(sprintf("  Root Mean Square Error: %.4f\n", x$goodness_of_fit$rmse))
  cat(sprintf("  Mean Absolute Percentage Error: %.2f%%\n", x$goodness_of_fit$mape * 100))
  cat(sprintf("  Correlation: %.4f\n", x$goodness_of_fit$correlation))
  cat("\n")
  
  # Statistical tests
  cat("STATISTICAL ADEQUACY TESTS:\n")
  
  # Ljung-Box test
  lb_test <- x$statistical_tests$ljung_box
  cat("  Ljung-Box Test for Residual Autocorrelation:\n")
  cat(sprintf("    X-squared = %.4f, df = %d, p-value = %.4f\n", 
              lb_test$statistic, lb_test$parameter, lb_test$p.value))
  cat("    H0: No residual autocorrelation (p > 0.05 indicates adequate model)\n")
  cat("    Result:", ifelse(lb_test$p.value > 0.05, " ADEQUATE", " INADEQUATE"), "\n\n")
  
  # Shapiro-Wilk test
  sw_test <- x$statistical_tests$shapiro_wilk
  cat("  Shapiro-Wilk Test for Normality:\n")
  cat(sprintf("    W = %.4f, p-value = %.4f\n", sw_test$statistic, sw_test$p.value))
  cat("    H0: Residuals are normally distributed (p > 0.05 indicates normality)\n")
  cat("    Result:", ifelse(sw_test$p.value > 0.05, " NORMAL", " NON-NORMAL"), "\n\n")
  
  # ARCH test
  arch_test <- x$statistical_tests$arch_effect
  cat("  ARCH Test for Heteroscedasticity:\n")
  cat(sprintf("    LM = %.4f, p-value = %.4f\n", arch_test$statistic, arch_test$p.value))
  cat("    H0: No ARCH effects (p > 0.05 indicates homoscedasticity)\n")
  cat("    Result:", ifelse(arch_test$p.value > 0.05, " HOMOSCEDASTIC", " HETEROSCEDASTIC"), "\n\n")
  
  # Kolmogorov-Smirnov test
  ks_test <- x$statistical_tests$distribution_similarity
  cat("  Kolmogorov-Smirnov Test for Distribution Similarity:\n")
  cat(sprintf("    D = %.4f, p-value = %.4f\n", ks_test$statistic, ks_test$p.value))
  cat("    H0: Original and synthetic data come from same distribution\n")
  cat("    Result:", ifelse(ks_test$p.value > 0.05, " SIMILAR DISTRIBUTIONS", " DIFFERENT DISTRIBUTIONS"), "\n\n")
  
  # Overall assessment
  cat("OVERALL ASSESSMENT:\n")
  adequate_tests <- sum(c(
    lb_test$p.value > 0.05,
    sw_test$p.value > 0.05,
    arch_test$p.value > 0.05,
    ks_test$p.value > 0.05
  ))
  cat(sprintf("  %d/4 tests passed adequacy criteria (p > 0.05)\n", adequate_tests))
  if (adequate_tests >= 3) {
    cat("  ✅ Model appears ADEQUATE for synthetic data generation\n")
  } else if (adequate_tests >= 2) {
    cat("  ⚠️  Model shows some INADEQUACIES\n")
  } else {
    cat("  ❌ Model appears INADEQUATE for synthetic data generation\n")
  }
}

# Helper functions for plotting
.plot_series_comparison <- function(x, which_sim, series_index, ...) {
  if (which_sim > x$n_sim) {
    stop("which_sim cannot exceed number of simulations")
  }
  
  if (x$is_multivariate) {
    if (series_index > ncol(x$original_series)) {
      stop("series_index cannot exceed number of series")
    }
    original <- x$original_series[, series_index]
    synthetic <- x$synthetic_series[[which_sim]][, series_index]
    series_name <- colnames(x$original_series)[series_index]
  } else {
    original <- x$original_series
    synthetic <- x$synthetic_series[[which_sim]]
    series_name <- "Series"
  }
  
  # Combine for common y-axis limits
  combined <- c(original, synthetic)
  ylim <- range(combined, na.rm = TRUE)
  
  par(mfrow = c(2, 1))
  
  # Plot original
  plot(original, main = paste("Original", series_name), 
       ylab = "Value", ylim = ylim, ...)
  grid()
  
  # Plot synthetic
  plot(synthetic, main = paste("Synthetic", series_name, "- Simulation", which_sim), 
       ylab = "Value", ylim = ylim, col = "red", ...)
  grid()
}

.plot_residuals <- function(x, series_index, ...) {
  if (x$is_multivariate) {
    if (series_index > ncol(x$original_series)) {
      stop("series_index cannot exceed number of series")
    }
    residuals <- na.omit(x$residuals[, series_index])
  } else {
    residuals <- na.omit(x$residuals)
  }
  
  par(mfrow = c(2, 2))
  
  # Time series of residuals
  plot(residuals, type = "l", main = "Residuals Time Series", 
       ylab = "Residuals", ...)
  abline(h = 0, col = "red", lty = 2)
  grid()
  
  # ACF of residuals
  acf(residuals, main = "ACF of Residuals", ...)
  
  # Histogram of residuals
  hist(residuals, main = "Distribution of Residuals", 
       xlab = "Residuals", probability = TRUE, ...)
  curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), 
        add = TRUE, col = "red", lwd = 2)
  
  # Q-Q plot
  qqnorm(residuals, main = "Q-Q Plot of Residuals", ...)
  qqline(residuals, col = "red")
}

.plot_acf_comparison <- function(x, which_sim, series_index, ...) {
  if (x$is_multivariate) {
    original <- x$original_series[, series_index]
    synthetic <- x$synthetic_series[[which_sim]][, series_index]
  } else {
    original <- x$original_series
    synthetic <- x$synthetic_series[[which_sim]]
  }
  
  par(mfrow = c(2, 1))
  acf(original, main = "ACF - Original Series", ...)
  acf(synthetic, main = "ACF - Synthetic Series", ...)
}

.plot_density_comparison <- function(x, which_sim, series_index, ...) {
  if (x$is_multivariate) {
    original <- na.omit(x$original_series[, series_index])
    synthetic <- na.omit(x$synthetic_series[[which_sim]][, series_index])
  } else {
    original <- na.omit(x$original_series)
    synthetic <- na.omit(x$synthetic_series[[which_sim]])
  }
  
  # Kernel density estimates
  dens_orig <- density(original)
  dens_synth <- density(synthetic)
  
  ylim <- range(c(dens_orig$y, dens_synth$y))
  xlim <- range(c(dens_orig$x, dens_synth$x))
  
  plot(dens_orig, main = "Density Comparison", xlim = xlim, ylim = ylim,
       xlab = "Value", col = "black", lwd = 2, ...)
  lines(dens_synth, col = "red", lwd = 2)
  legend("topright", legend = c("Original", "Synthetic"), 
         col = c("black", "red"), lwd = 2)
  grid()
}

.plot_qq_comparison <- function(x, which_sim, series_index, ...) {
  if (x$is_multivariate) {
    original <- na.omit(x$original_series[, series_index])
    synthetic <- na.omit(x$synthetic_series[[which_sim]][, series_index])
  } else {
    original <- na.omit(x$original_series)
    synthetic <- na.omit(x$synthetic_series[[which_sim]])
  }
  
  qqplot(original, synthetic, main = "Q-Q Plot: Original vs Synthetic",
         xlab = "Original Quantiles", ylab = "Synthetic Quantiles", ...)
  abline(0, 1, col = "red", lwd = 2)
  grid()
}

# Helper functions for summary
.get_basic_info <- function(x) {
  list(
    original_length = if (x$is_multivariate) nrow(x$original_series) else length(x$original_series),
    synthetic_length = x$length_out,
    n_sim = x$n_sim,
    frequency = x$ts_properties$frequency,
    is_multivariate = x$is_multivariate,
    n_series = if (x$is_multivariate) ncol(x$original_series) else 1
  )
}

.get_model_info <- function(x) {
  list(
    model_type = x$model_info$model_type,
    formula = x$model_info$formula,
    order = x$model_info$order
  )
}

.get_goodness_of_fit <- function(x) {
  # Compare first synthetic series with original
  if (x$is_multivariate) {
    # Use first series for comparison
    original <- x$original_series[, 1]
    synthetic <- x$synthetic_series[[1]][, 1]
  } else {
    original <- x$original_series
    synthetic <- x$synthetic_series[[1]]
  }
  
  # Ensure same length for comparison
  min_len <- min(length(original), length(synthetic))
  original_trim <- original[1:min_len]
  synthetic_trim <- synthetic[1:min_len]
  
  mae <- mean(abs(original_trim - synthetic_trim), na.rm = TRUE)
  rmse <- sqrt(mean((original_trim - synthetic_trim)^2, na.rm = TRUE))
  mape <- mean(abs((original_trim - synthetic_trim) / original_trim), na.rm = TRUE)
  correlation <- cor(original_trim, synthetic_trim, use = "complete.obs")
  
  list(mae = mae, rmse = rmse, mape = mape, correlation = correlation)
}

.perform_adequacy_tests <- function(x) {
  # Use residuals from the first series for multivariate case
  if (x$is_multivariate) {
    residuals <- na.omit(x$residuals[, 1])
    original <- na.omit(x$original_series[, 1])
    synthetic <- na.omit(x$synthetic_series[[1]][, 1])
  } else {
    residuals <- na.omit(x$residuals)
    original <- na.omit(x$original_series)
    synthetic <- na.omit(x$synthetic_series[[1]])
  }
  
  # 1. Ljung-Box test for residual autocorrelation
  lb_test <- Box.test(residuals, lag = min(20, length(residuals) %/% 4), type = "Ljung-Box")
  
  # 2. Shapiro-Wilk test for normality
  sw_test <- shapiro.test(residuals)
  
  # 3. ARCH test for heteroscedasticity (using LM test)
  arch_test <- .arch_lm_test(residuals)
  
  # 4. Kolmogorov-Smirnov test for distribution similarity
  ks_test <- ks.test(original, synthetic)
  
  list(
    ljung_box = lb_test,
    shapiro_wilk = sw_test,
    arch_effect = arch_test,
    distribution_similarity = ks_test
  )
}

.arch_lm_test <- function(residuals, lag = 5) {
  # Lagrange Multiplier test for ARCH effects
  n <- length(residuals)
  squared_resid <- residuals^2
  
  # Create lag matrix
  lag_matrix <- sapply(1:lag, function(i) c(rep(NA, i), squared_resid[1:(n - i)]))
  
  # Remove NA rows
  complete_cases <- complete.cases(lag_matrix)
  y <- squared_resid[(lag + 1):n]
  X <- cbind(1, lag_matrix[(lag + 1):n, ])
  
  # Fit linear model
  fit <- lm(y ~ X - 1)
  r_squared <- summary(fit)$r.squared
  
  # LM statistic
  lm_stat <- n * r_squared
  p_value <- pchisq(lm_stat, df = lag, lower.tail = FALSE)
  
  list(statistic = lm_stat, p.value = p_value, method = "ARCH LM Test")
}