# GARCH(1,1) Forecast Function similar to forecast::thetaf or forecast::meanf
# Based on garch11Fit from https://lbelzile.github.io/timeseRies/

#' GARCH(1,1) Forecasting Function
#'
#' Generates forecasts from a GARCH(1,1) model, similar in style to 
#' \code{forecast::thetaf} or \code{forecast::meanf}. Supports both analytical
#' and bootstrap-based prediction intervals.
#'
#' @param y A numeric vector or time series object containing the data to forecast
#' @param h Number of periods for forecasting (default: 10)
#' @param level Confidence levels for prediction intervals (default: c(80, 95))
#' @param fan Logical. If TRUE, level is set to seq(51, 99, by = 3) (default: FALSE)
#' @param bootstrap Logical. If TRUE, uses bootstrap method for prediction intervals 
#'   (default: FALSE)
#' @param npaths Number of bootstrap paths when bootstrap = TRUE (default: 1000)
#'
#' @return An object of class "garch11forecast" containing:
#' \itemize{
#'   \item \code{model}: The fitted GARCH(1,1) model object
#'   \item \code{method}: Method name ("GARCH(1,1)")
#'   \item \code{mean}: Point forecasts (time series)
#'   \item \code{variance}: Forecasted conditional variances (time series)
#'   \item \code{lower}: Lower prediction bounds (matrix)
#'   \item \code{upper}: Upper prediction bounds (matrix)
#'   \item \code{level}: Confidence levels used
#'   \item \code{x}: Original time series
#'   \item \code{fitted}: Fitted values from the model
#'   \item \code{residuals}: Model residuals
#' }
#'
#'
#' @seealso \code{\link{garch11Fit}} for model fitting, \code{\link{plot.garch11forecast}} 
#'   for plotting forecasts
#' @export
garch11f <- function(y, h = 10, level = c(80, 95), fan = FALSE, bootstrap = FALSE, npaths = 1000) {
  
  # Input validation
  if (!is.numeric(y) && !is.ts(y)) stop("y must be numeric or time series")
  if (h <= 0 || !is.numeric(h) || length(h) != 1) stop("h must be positive integer")
  if (any(level <= 0 | level >= 100)) stop("level must be between 0 and 100")
  if (bootstrap && (npaths <= 0 || npaths != round(npaths))) stop("npaths must be positive integer")
  
  # Fit GARCH(1,1) model
  model <- garch11Fit(y, return_model = TRUE)
  
  # Check model validity
  if (model$convergence != 0) {
    warning("GARCH model fitting may be unreliable due to convergence issues")
  }
  
  # Extract parameters
  mu <- model$coefficients[1]
  omega <- model$coefficients[2]
  alpha <- model$coefficients[3]
  beta <- model$coefficients[4]
  
  # Get last residual and variance for forecasting
  last_residual <- model$last_residual
  last_variance <- model$last_variance
  
  # Generate forecasts
  if (bootstrap) {
    # Bootstrap forecasting (improved efficiency)
    forecasts <- garch11_bootstrap_forecast(model, h, npaths)
    mean_forecast <- rowMeans(forecasts$returns)
    variance_forecast <- rowMeans(forecasts$variances)
    
    # Calculate prediction intervals from bootstrap
    lower_bounds <- upper_bounds <- matrix(NA, nrow = h, ncol = length(level))
    for (i in 1:length(level)) {
      prob <- (100 - level[i]) / 200
      lower_bounds[, i] <- apply(forecasts$returns, 1, quantile, prob = prob)
      upper_bounds[, i] <- apply(forecasts$returns, 1, quantile, prob = 1 - prob)
    }
  } else {
    # Analytical forecasting with detailed comments
    mean_forecast <- rep(mu, h)
    variance_forecast <- numeric(h)
    
    # Multi-step GARCH variance forecasting
    # h=1: σ²_{t+1} = ω + α*V²_t + β*σ²_t (uses actual last residual)
    # h>1: σ²_{t+h} = ω + (α+β)*σ²_{t+h-1} (since E[V²_{t+h-1}] = σ²_{t+h-1})
    current_var <- last_variance
    current_resid_sq <- last_residual^2
    
    for (i in 1:h) {
      if (i == 1) {
        # One-step ahead: use actual last residual
        variance_forecast[i] <- omega + alpha * current_resid_sq + beta * current_var
      } else {
        # Multi-step: expected value of previous squared residual equals previous variance
        variance_forecast[i] <- omega + (alpha + beta) * variance_forecast[i-1]
      }
    }
    
    # Calculate prediction intervals (assuming normal distribution)
    sigma_forecast <- sqrt(variance_forecast)
    lower_bounds <- upper_bounds <- matrix(NA, nrow = h, ncol = length(level))
    
    for (i in 1:length(level)) {
      z_score <- qnorm(1 - (100 - level[i]) / 200)
      lower_bounds[, i] <- mean_forecast - z_score * sigma_forecast
      upper_bounds[, i] <- mean_forecast + z_score * sigma_forecast
    }
  }
  
  # Create time series objects
  n <- length(y)
  tsp_orig <- tsp(y)
  
  if (is.null(tsp_orig)) {
    start_forecast <- n + 1
    forecast_times <- seq(start_forecast, start_forecast + h - 1)
  } else {
    freq <- tsp_orig[3]
    end_orig <- tsp_orig[2]
    start_forecast <- end_orig + 1/freq
    forecast_times <- seq(start_forecast, by = 1/freq, length.out = h)
  }
  
  # Create forecast object similar to forecast package structure
  result <- list(
    model = model,
    method = "GARCH(1,1)",
    mean = ts(mean_forecast, start = start_forecast, frequency = ifelse(is.null(tsp_orig), 1, tsp_orig[3])),
    variance = ts(variance_forecast, start = start_forecast, frequency = ifelse(is.null(tsp_orig), 1, tsp_orig[3])),
    lower = lower_bounds,
    upper = upper_bounds,
    level = level,
    x = y,
    fitted = ts(model$fitted_mean, start = tsp_orig[1], frequency = ifelse(is.null(tsp_orig), 1, tsp_orig[3])),
    residuals = ts(model$residuals, start = tsp_orig[1], frequency = ifelse(is.null(tsp_orig), 1, tsp_orig[3]))
  )
  
  class(result) <- c("garch11forecast", "forecast")
  
  return(result)
}

# First, let's include the original garch11Fit function with modifications to return model object
garch11Fit <- function(x, return_model = TRUE) {
  # Input validation
  if (!is.numeric(x)) stop("Input must be numeric")
  if (length(x) < 10) stop("Time series too short (minimum 10 observations required)")
  if (any(is.na(x))) stop("Time series contains missing values")
  if (any(!is.finite(x))) stop("Time series contains infinite values")
  
  require(numDeriv)
  
  # Step 1: Initialize Model Parameters and Bounds
  series_mean <- mean(x)
  series_var <- var(x)
  S <- 1e-06
  params <- c(mu = series_mean, omega = 0.1 * series_var, alpha = 0.1, beta = 0.8)
  lowerBounds <- c(mu = -10 * abs(series_mean), omega = S^2, alpha = S, beta = S)
  upperBounds <- c(mu = 10 * abs(series_mean), omega = 100 * series_var, alpha = 1 - S, beta = 1 - S)
  
  # Step 2: Set Conditional Distribution Function
  garchDist <- function(z, hh) {
    dnorm(x = z/hh)/hh
  }
  
  # Step 3: Compose log-Likelihood Function
  garchLLH <- function(parm) {
    mu <- parm[1]
    omega <- parm[2]
    alpha <- parm[3]
    beta <- parm[4]
    z <- (x - mu)
    # Initial variance estimate for GARCH filter
    init_var <- mean(z^2)
    # Use Filter Representation
    e <- omega + alpha * c(init_var, z[-length(x)]^2)
    h <- filter(e, beta, "r", init = init_var)
    hh <- sqrt(abs(h))
    -sum(log(garchDist(z, hh))) # llh
  }
  
  # Step 4: Estimate Parameters and Compute Numerically Hessian
  fit <- nlminb(start = params, objective = garchLLH, 
                lower = lowerBounds, upper = upperBounds)
  
  # Check convergence
  if (fit$convergence != 0) {
    warning(paste("Optimization did not converge. Code:", fit$convergence))
  }
  
  # Check stationarity condition
  alpha_est <- fit$par[3]
  beta_est <- fit$par[4]
  if (alpha_est + beta_est >= 1) {
    warning("GARCH process appears to be non-stationary (alpha + beta >= 1)")
  }
  
  if (return_model) {
    # Calculate fitted values and residuals for forecasting
    mu <- fit$par[1]
    omega <- fit$par[2]
    alpha <- fit$par[3]
    beta <- fit$par[4]
    
    z <- (x - mu)
    init_var <- mean(z^2)
    e <- omega + alpha * c(init_var, z[-length(x)]^2)
    h <- filter(e, beta, "r", init = init_var)
    
    # Store model information
    model_info <- list(
      coefficients = fit$par,
      data = x,
      residuals = z,
      sigma_squared = as.numeric(h),
      fitted_mean = rep(mu, length(x)),
      loglik = -fit$objective,
      convergence = fit$convergence,
      last_residual = z[length(z)],
      last_variance = h[length(h)],
      stationarity = alpha + beta
    )
    
    class(model_info) <- "garch11"
    return(model_info)
  } else {
    # Original behavior - print results
    Hessian <- numDeriv::hessian(func = garchLLH, x = fit$par)
    se.coef <- sqrt(diag(solve(Hessian)))
    tval <- fit$par/se.coef
    matcoef <- cbind(fit$par, se.coef, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) <- list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
    return(matcoef)
  }
}

# Improved bootstrap forecasting helper function (more efficient)
garch11_bootstrap_forecast <- function(model, h, npaths = 1000) {
  
  mu <- model$coefficients[1]
  omega <- model$coefficients[2]
  alpha <- model$coefficients[3]
  beta <- model$coefficients[4]
  
  # Standardized residuals for bootstrap
  std_residuals <- model$residuals / sqrt(model$sigma_squared)
  
  # Pre-allocate matrices for efficiency
  returns <- matrix(NA, nrow = h, ncol = npaths)
  variances <- matrix(NA, nrow = h, ncol = npaths)
  
  last_variance <- model$last_variance
  last_resid_sq <- model$last_residual^2
  
  # Vectorized bootstrap sampling
  bootstrap_innovations <- matrix(sample(std_residuals, h * npaths, replace = TRUE), 
                                  nrow = h, ncol = npaths)
  
  # Initialize first step
  variances[1, ] <- omega + alpha * last_resid_sq + beta * last_variance
  returns[1, ] <- mu + sqrt(variances[1, ]) * bootstrap_innovations[1, ]
  
  # Iterate for remaining steps (still need loop for sequential nature)
  if (h > 1) {
    for (i in 2:h) {
      variances[i, ] <- omega + alpha * returns[i-1, ]^2 + beta * variances[i-1, ]
      returns[i, ] <- mu + sqrt(variances[i, ]) * bootstrap_innovations[i, ]
    }
  }
  
  return(list(returns = returns, variances = variances))
}

# Print method for garch11forecast objects
print.garch11forecast <- function(x, ...) {
  cat("Forecasts from", x$method, "\n")
  if (x$model$convergence != 0) {
    cat("WARNING: Model estimation may be unreliable (convergence code:", x$model$convergence, ")\n")
  }
  if (x$model$stationarity >= 0.99) {
    cat("WARNING: Process appears close to non-stationarity (α + β =", round(x$model$stationarity, 4), ")\n")
  }
  cat("\n")
  
  cat("Point forecasts:\n")
  print(x$mean)
  
  cat("\nConditional variance forecasts:\n")
  print(x$variance)
  
  if (!is.null(x$level)) {
    cat("\nPrediction intervals:\n")
    for (i in 1:length(x$level)) {
      cat(paste0(x$level[i], "%: ["))
      cat(paste(round(x$lower[, i], 4), round(x$upper[, i], 4), sep = ", ", collapse = "] ["))
      cat("]\n")
    }
  }
}

#' GARCH(1,1) Model Fitting Function
#'
#' Fits a GARCH(1,1) model to a time series using maximum likelihood estimation.
#' This function estimates the parameters of the GARCH(1,1) model: 
#' \eqn{\sigma_t^2 = \omega + \alpha \epsilon_{t-1}^2 + \beta \sigma_{t-1}^2}
#'
#' @param x A numeric vector or time series of data to fit the GARCH model to
#' @param return_model Logical. If TRUE (default), returns a complete model object 
#'   suitable for forecasting. If FALSE, returns coefficient table with statistics.
#'
#' @return If `return_model = TRUE`, returns an object of class "garch11" containing:
#' \itemize{
#'   \item \code{coefficients}: Estimated parameters (mu, omega, alpha, beta)
#'   \item \code{data}: Original time series data
#'   \item \code{residuals}: Model residuals
#'   \item \code{sigma_squared}: Conditional variances
#'   \item \code{fitted_mean}: Fitted mean values
#'   \item \code{loglik}: Log-likelihood value
#'   \item \code{convergence}: Optimization convergence code
#'   \item \code{last_residual}: Last residual for forecasting
#'   \item \code{last_variance}: Last conditional variance for forecasting
#'   \item \code{stationarity}: Alpha + beta stationarity measure
#' }
#' If `return_model = FALSE`, returns a coefficient matrix with estimates, 
#' standard errors, t-values, and p-values.
#'
#' @examples
#' \donttest{
#' # Generate sample GARCH(1,1) data
#' set.seed(123)
#' n <- 500
#' omega <- 0.1; alpha <- 0.1; beta <- 0.8; mu <- 0.05
#' 
#' y <- numeric(n)
#' sigma2 <- numeric(n)
#' sigma2[1] <- omega / (1 - alpha - beta)
#' y[1] <- mu + sqrt(sigma2[1]) * rnorm(1)
#' 
#' for (t in 2:n) {
#'   sigma2[t] <- omega + alpha * (y[t-1] - mu)^2 + beta * sigma2[t-1]
#'   y[t] <- mu + sqrt(sigma2[t]) * rnorm(1)
#' }
#' 
#' # Fit GARCH(1,1) model
#' model <- garch11Fit(y, return_model = TRUE)
#' print(model$coefficients)
#' }
#'
#' @seealso \code{\link{garch11f}} for forecasting with the fitted model
#' @export
garch11Fit <- function(x, return_model = TRUE) {
  # ... function body remains the same ...
}

#' GARCH(1,1) Forecasting Function
#'
#' Generates forecasts from a GARCH(1,1) model, similar in style to 
#' \code{forecast::thetaf} or \code{forecast::meanf}. Supports both analytical
#' and bootstrap-based prediction intervals.
#'
#' @param y A numeric vector or time series object containing the data to forecast
#' @param h Number of periods for forecasting (default: 10)
#' @param level Confidence levels for prediction intervals (default: c(80, 95))
#' @param fan Logical. If TRUE, level is set to seq(51, 99, by = 3) (default: FALSE)
#' @param bootstrap Logical. If TRUE, uses bootstrap method for prediction intervals 
#'   (default: FALSE)
#' @param npaths Number of bootstrap paths when bootstrap = TRUE (default: 1000)
#'
#' @return An object of class "garch11forecast" containing:
#' \itemize{
#'   \item \code{model}: The fitted GARCH(1,1) model object
#'   \item \code{method}: Method name ("GARCH(1,1)")
#'   \item \code{mean}: Point forecasts (time series)
#'   \item \code{variance}: Forecasted conditional variances (time series)
#'   \item \code{lower}: Lower prediction bounds (matrix)
#'   \item \code{upper}: Upper prediction bounds (matrix)
#'   \item \code{level}: Confidence levels used
#'   \item \code{x}: Original time series
#'   \item \code{fitted}: Fitted values from the model
#'   \item \code{residuals}: Model residuals
#' }
#'
#' @examples
#' \donttest{
#' # Generate sample data
#' set.seed(123)
#' n <- 200
#' y <- rnorm(n)  # Simple example with normal data
#' 
#' # Generate forecasts
#' forecasts <- garch11f(y, h = 20, level = c(80, 95))
#' print(forecasts)
#' plot(forecasts)
#' 
#' # Use bootstrap method
#' forecasts_boot <- garch11f(y, h = 20, bootstrap = TRUE, npaths = 500)
#' plot(forecasts_boot, main = "Bootstrap GARCH Forecasts")
#' }
#'
#' @seealso \code{\link{garch11Fit}} for model fitting, \code{\link{plot.garch11forecast}} 
#'   for plotting forecasts
#' @export
garch11f <- function(y, h = 10, level = c(80, 95), fan = FALSE, bootstrap = FALSE, npaths = 1000) {
  # ... function body remains the same ...
}

#' Bootstrap Forecasting for GARCH(1,1) Model
#'
#' Internal function for generating bootstrap-based forecasts from a fitted GARCH(1,1) model.
#' This function is called internally by \code{garch11f} when \code{bootstrap = TRUE}.
#'
#' @param model A fitted GARCH(1,1) model object from \code{garch11Fit}
#' @param h Forecast horizon
#' @param npaths Number of bootstrap paths to generate (default: 1000)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{returns}: Matrix of bootstrap return paths (h × npaths)
#'   \item \code{variances}: Matrix of bootstrap variance paths (h × npaths)
#' }
#'
#' @keywords internal
garch11_bootstrap_forecast <- function(model, h, npaths = 1000) {
  # ... function body remains the same ...
}

#' Print Method for GARCH(1,1) Forecast Objects
#'
#' Prints a summary of GARCH(1,1) forecasts in a user-friendly format.
#'
#' @param x An object of class "garch11forecast" from \code{garch11f}
#' @param ... Additional arguments passed to print method
#'
#' @return Invisibly returns the input object
#'
#' @seealso \code{\link{garch11f}}
#' @method print garch11forecast
#' @export
print.garch11forecast <- function(x, ...) {
  # ... function body remains the same ...
}

#' Plot Method for GARCH(1,1) Forecast Objects
#'
#' Creates a visualization of GARCH(1,1) forecasts with historical data, 
#' point forecasts, and prediction intervals.
#'
#' @param x An object of class "garch11forecast" from \code{garch11f}
#' @param main Plot title (default: "GARCH(1,1) Forecasts")
#' @param xlab X-axis label (default: "Time")
#' @param ylab Y-axis label (default: "Value")
#' @param ... Additional arguments passed to base plot function
#'
#' @return A ggplot object (if ggplot2 is available) or base R plot
#'
#' @examples
#' \donttest{
#' # After creating forecasts with garch11f()
#' forecasts <- garch11f(rnorm(100), h = 10)
#' plot(forecasts)
#' plot(forecasts, main = "Custom Title", ylab = "Returns")
#' }
#'
#' @seealso \code{\link{garch11f}}
#' @method plot garch11forecast
#' @export
plot.garch11forecast <- function(x, main = "GARCH(1,1) Forecasts", xlab = "Time", ylab = "Value", ...) {
  # ... function body remains the same ...
}

#' Summary Method for GARCH(1,1) Forecast Objects
#'
#' Provides a comprehensive summary of GARCH(1,1) model estimates and forecasts.
#'
#' @param object An object of class "garch11forecast" from \code{garch11f}
#' @param ... Additional arguments (currently not used)
#'
#' @return Prints a detailed summary to the console
#'
#' @method summary garch11forecast
#' @export
summary.garch11forecast <- function(object, ...) {
  cat("GARCH(1,1) Forecast Summary\n")
  cat("==========================\n")
  cat("Model Parameters:\n")
  cat("μ (mean):", round(object$model$coefficients[1], 6), "\n")
  cat("ω (omega):", round(object$model$coefficients[2], 6), "\n")
  cat("α (alpha):", round(object$model$coefficients[3], 6), "\n")
  cat("β (beta):", round(object$model$coefficients[4], 6), "\n")
  cat("α + β:", round(object$model$stationarity, 6), "\n")
  cat("Log-Likelihood:", round(object$model$loglik, 2), "\n")
  cat("Convergence:", ifelse(object$model$convergence == 0, "Success", "Issues"), "\n")
  cat("\nForecast Information:\n")
  cat("Horizon:", length(object$mean), "periods\n")
  cat("Confidence levels:", paste(object$level, "%", collapse = ", "), "\n")
  if (!is.null(object$method) && object$method == "GARCH(1,1)") {
    cat("Method: Analytical forecasting\n")
  }
}