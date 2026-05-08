#' Model-agnostic statistical probabilistic forecasting with conditional volatility
#'
#' Combines a mean model and a volatility model to produce probabilistic
#' forecasts. Uncertainty is propagated from two sources: simulated innovation
#' draws and simulated volatility paths derived from the variance model's
#' residuals via bootstrap.
#'
#' @param y A numeric vector or time series of class \code{ts}.
#' @param h Forecasting horizon.
#' @param mean_model Function to fit the mean model (default: \code{forecast::auto.arima}).
#'   Options include \code{forecast::auto.arima}, \code{forecast::ets}, \code{forecast::thetaf},
#'   or any custom function that returns an object with \code{forecast} method.
#' @param sigma_model Function to fit the variance model on squared residuals
#'   (default: \code{forecast::auto.arima}).
#' @param innovation Distribution for standardized innovations. One of
#'   \code{"gaussian"} (fastest, often sufficient), \code{"student"} (heavy-tailed),
#'   or \code{"empirical"} (non-parametric, most flexible).
#' @param level Confidence level for prediction intervals (default: 95).
#' @param B Number of simulated paths (default: 2000).
#' @param bootstrap_vol Whether to bootstrap volatility residuals (default: TRUE).
#'   If FALSE, uses parametric normal approximation for volatility.
#' @param ... Additional arguments passed to \code{mean_model} and \code{sigma_model}.
#'
#' @return A \code{forecast} object with components:
#'   \describe{
#'     \item{\code{x}}{Original time series.}
#'     \item{\code{mean}}{Point forecast (mean of simulations).}
#'     \item{\code{lower}}{Lower prediction interval bound.}
#'     \item{\code{upper}}{Upper prediction interval bound.}
#'     \item{\code{level}}{Confidence level.}
#'     \item{\code{sims}}{Matrix of simulated forecast paths (\code{h x B}).}
#'     \item{\code{model}}{List containing \code{mean} and \code{sigma} model fits.}
#'     \item{\code{residuals}}{Residuals from the mean model.}
#'     \item{\code{standardized_residuals}}{Scaled standardized residuals used for innovation fitting.}
#'     \item{\code{method}}{Character string describing the method used.}
#'   }
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(forecast)
#'
#' # Basic usage with Google stock data
#' y <- fpp2::goog200
#'
#' # Gaussian innovations (fastest, often sufficient)
#' fc1 <- condvolf(y, h = 20, innovation = "gaussian")
#' plot(fc1)
#'
#' # Student-t innovations (heavy-tailed)
#' fc2 <- condvolf(y, h = 20, innovation = "student")
#' plot(fc2)
#'
#' # Empirical innovations (non-parametric)
#' fc3 <- condvolf(y, h = 20, innovation = "empirical")
#' plot(fc3)
#'
#' # Using different mean and volatility models
#' fc4 <- condvolf(y, h = 20,
#'                 mean_model = forecast::thetaf,
#'                 sigma_model = forecast::ets,
#'                 innovation = "gaussian")
#'
#' # Compare prediction intervals
#' par(mfrow = c(2, 2))
#' plot(fc1, main = "Gaussian innovations")
#' plot(fc2, main = "Student-t innovations")
#' plot(fc3, main = "Empirical innovations")
#' plot(fc4, main = "Different mean and volatility models")
#' par(mfrow = c(1, 1))
#' }
condvolf <- function(
    y,
    h = 10L,
    mean_model = forecast::auto.arima,
    sigma_model = forecast::auto.arima,
    innovation = c("empirical", "gaussian", "student"),
    level = 95,
    B = 1000L,
    bootstrap_vol = TRUE,
    ...
)
{
  # Match innovation type
  innovation <- match.arg(innovation)
  
  # Ensure y is a time series
  if (!is.ts(y)) {
    y <- ts(y)
  }
  
  # Store original data and metadata
  x <- y
  freq_x <- frequency(y)
  start_preds <- tsp(y)[2] + 1 / freq_x
  
  #---------------------------------------------------------------------------
  # 1. Fit mean model and obtain forecasts
  #---------------------------------------------------------------------------
  
  # Use generic forecast function that works with auto.arima, ets, thetaf, etc.
  fit_mean <- ahead::genericforecast(
    FUN = mean_model,
    y = y,
    h = h,
    ...
  )
  
  # Extract point forecasts and residuals
  mu_hat <- as.numeric(fit_mean$mean)
  resids <- as.numeric(residuals(fit_mean))
  
  # Handle potential NA residuals
  if (any(is.na(resids))) {
    resids <- na.omit(resids)
    if (length(resids) == 0) {
      resids <- rep(0, length(y))
    }
  }
  
  #---------------------------------------------------------------------------
  # 2. Fit volatility model on squared residuals
  #---------------------------------------------------------------------------
  
  fit_sigma <- ahead::genericforecast(
    FUN = sigma_model,
    y = resids^2,
    h = h,
    ...
  )
  
  # In-sample fitted volatility
  sigma_insample <- sqrt(pmax(as.numeric(fitted(fit_sigma)), 1e-8))
  
  # Standardized residuals (innovation sample)
  z <- as.numeric(scale(resids / sigma_insample))
  z[is.na(z) | !is.finite(z)] <- 0
  
  #---------------------------------------------------------------------------
  # 3. Simulate volatility paths
  #---------------------------------------------------------------------------
  
  var_mean <- pmax(as.numeric(fit_sigma$mean), 1e-8)
  
  if (bootstrap_vol) {
    # Bootstrap from volatility model residuals (recommended)
    var_residuals <- as.numeric(residuals(fit_sigma))
    var_residuals <- var_residuals[is.finite(var_residuals)]
    
    if (length(var_residuals) == 0) {
      var_residuals <- 0
    }
    
    var_sims <- matrix(
      sample(var_residuals, size = h * B, replace = TRUE),
      nrow = h, ncol = B
    )
  } else { # if !bootstrap_vol
    # Parametric simulation (fallback)
    var_sims <- matrix(0, nrow = h, ncol = B)
  }
  
  sigma_sims <- sqrt(pmax(var_mean + var_sims, 1e-8))
  
  #---------------------------------------------------------------------------
  # 4. Simulate innovations
  #---------------------------------------------------------------------------
  
  zsim <- matrix(0, nrow = h, ncol = B)
  
  if (innovation == "gaussian") {
    
    # Gaussian (normal) innovations
    zsim <- matrix(rnorm(h * B), nrow = h, ncol = B)
    
  } else if (innovation == "student") {
    
    # Student-t innovations with fitted degrees of freedom
    fit_t <- tryCatch(
      MASS::fitdistr(z, densfun = "t"),
      error = function(e) NULL
    )
    
    if (is.null(fit_t)) {
      warning("Student-t fit failed; falling back to Gaussian innovations.")
      zsim <- matrix(rnorm(h * B), nrow = h, ncol = B)
    } else {
      df_t <- max(fit_t$estimate["df"], 3)  # Ensure df >= 3 for finite variance
      zsim <- matrix(rt(h * B, df = df_t), nrow = h, ncol = B)
    }
    
  } else { # innovation == "empirical"
    # Non-parametric bootstrap from standardized residuals
    if (length(z) == 0) {
      warning("No standardized residuals available; using Gaussian innovations.")
      zsim <- matrix(rnorm(h * B), nrow = h, ncol = B)
    } else {
      zsim <- matrix(
        sample(z, size = h * B, replace = TRUE),
        nrow = h, ncol = B
      )
    }
  }
  
  #---------------------------------------------------------------------------
  # 5. Reconstruct forecast paths
  #---------------------------------------------------------------------------
  
  # Combine mean, volatility, and innovations
  f <- matrix(mu_hat, nrow = h, ncol = B) + zsim * sigma_sims
  
  #---------------------------------------------------------------------------
  # 6. Summarize forecasts
  #---------------------------------------------------------------------------
  
  # Point forecast: mean of simulations
  mean_f <- rowMeans(f)
  
  # Prediction intervals: quantiles of simulations
  alpha <- 1 - level / 100
  lower <- apply(f, 1, quantile, probs = alpha / 2, na.rm = TRUE)
  upper <- apply(f, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  
  #---------------------------------------------------------------------------
  # 7. Prepare output
  #---------------------------------------------------------------------------
  
  out <- list()
  out$x <- x
  out$mean <- ts(mean_f, start = start_preds, frequency = freq_x)
  out$lower <- ts(lower, start = start_preds, frequency = freq_x)
  out$upper <- ts(upper, start = start_preds, frequency = freq_x)
  out$level <- level
  out$sims <- f  # Matrix of simulations
  out$model <- list(mean = fit_mean, sigma = fit_sigma)
  out$residuals <- resids
  out$standardized_residuals <- z
  out$method <- sprintf("condvolf (mean=%s, sigma=%s, innov=%s)",
                        deparse(substitute(mean_model)),
                        deparse(substitute(sigma_model)),
                        innovation)
  
  # Ensure proper class for forecast package compatibility
  class(out) <- c("condvolf", "forecast")
  
  return(out)
}
