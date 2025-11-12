#' Context-Aware Theta method forecast
#'
#' Returns forecasts and prediction intervals for a context-aware theta method forecast.
#'
#' The classical theta method of Assimakopoulos and Nikolopoulos (2000) is equivalent to
#' simple exponential smoothing with drift, using theta = 2. This implementation extends
#' the method by:
#' 1. Allowing flexible specification of theta parameter (default 0.5, which gives theta = 2)
#' 2. Using machine learning models to capture non-linear drift patterns
#' 3. Computing time-varying slopes via numerical differentiation
#'
#' @param y Time series (ts object or numeric vector)
#' @param h Forecast horizon (default: 2*frequency for seasonal, 10 for non-seasonal)
#' @param level Confidence level(s) for prediction intervals (default: 95)
#' @param theta Theta coefficient for drift weighting (default: 0.5, giving theta line = 2)
#'   - theta = 0: No drift (pure SES)
#'   - theta = 0.5: Classical Theta method (theta line = 2)
#'   - theta = 1: Full drift
#'   - theta > 1: Amplified drift
#' @param fit_func Model fitting function (default: lm)
#' @param predict_func Prediction function (default: predict)
#' @param x Deprecated, use y instead
#' @param ... Additional arguments passed to fit_func
#' @param type_pi 
#' @param nsim 
#' @param block_size 
#' @param seed 
#' @param B 
#'
#' @return An object of class `forecast`
#' @author T. Moudiki (based on RJH)
#' @seealso [forecast::thetaf()], [forecast::ses()]
#' @references
#'
#' Hyndman, R.J., and Billah, B. (2003) Unmasking the Theta method.
#' \emph{International J. Forecasting}, \bold{19}, 287-290.
#'
#' @keywords ts
#' @examples
#' # Classical theta-like behavior (theta = 0.5)
#' fit1 <- ctxthetaf(AirPassengers)
#' plot(fit1)
#'
#' # No drift (theta = 0)
#' fit2 <- ctxthetaf(AirPassengers, theta = 0)
#' plot(fit2)
#'
#' # Amplified drift (theta = 1)
#' fit3 <- ctxthetaf(AirPassengers, theta = 1)
#' plot(fit3)
#'
#' # With Random Forest
#' fit4 <- ctxthetaf(AirPassengers, fit_func = randomForest::randomForest)
#' plot(fit4)
#'
#' @export
ctxthetaf <- function(y,
                      h = if (frequency(y) > 1)
                        2 * frequency(y)
                      else
                        10,
                      level = 95,
                      theta = 0.5,
                      fit_func = lm,
                      predict_func = predict,
                      type_pi = c("gaussian", "block-bootstrap", 
                                 "surrogate", "kde", "bootstrap", 
                                 "fitdistr", "meboot"),
                      nsim = 100L, 
                      block_size = 5,
                      seed = 123L, 
                      B = NULL,
                      x = y,
                      ...) {
  stopifnot(h > 0)
  stopifnot(is.numeric(theta))
  type_pi <- match.arg(type_pi)
  
  if (type_pi != "gaussian")
  {
    fcast_func <- function(y, h, level)
    {
      return(ahead::ctxthetaf(y = y,
                              h = h,
                              level = level,
                              theta = theta,
                              fit_func = fit_func,
                              predict_func = predict_func,
                              type_pi = "gaussian",
                              x = y,
                              ...))
    }
    return(ahead::conformalize(fcast_func, 
                 y = y, 
                 h = h, 
                 level=level,
                 method = type_pi,
                 nsim = nsim, 
                 block_size = block_size,
                 seed = seed, 
                 B = B,
                 ...))
  }
  
  # Check seasonality
  n <- length(x)
  x <- as.ts(x)
  m <- frequency(x)
  if (m > 1 && !is.constant(x) && n > 2 * m) {
    r <- as.numeric(acf(x, lag.max = m, plot = FALSE)$acf)[-1]
    stat <- sqrt((1 + 2 * sum(r[-m]^2)) / n)
    seasonal <- (abs(r[m]) / stat > qnorm(0.95))
  } else {
    seasonal <- FALSE
  }
  
  # Seasonal decomposition
  origx <- x
  if (seasonal) {
    decomp <- decompose(x, type = "multiplicative")
    if (any(abs(seasonal(decomp)) < 1e-4)) {
      warning("Seasonal indexes close to zero. Using non-seasonal Theta method")
    } else {
      x <- seasadj(decomp)
    }
  }
  
  # Find theta lines
  fcast <- forecast::ses(x, h = h)
  
  # Get slopes for historical + forecast horizon
  all_slopes <- theta * estimate_theta_slope(fit_func, predict_func, y, h = h, ...)
  
  # Use the LAST h slopes (the future ones)
  tmp2 <- tail(all_slopes, h)
  
  # Apply drift adjustment
  alpha <- pmax(1e-10, fcast$model$par["alpha"])
  fcast$mean <- fcast$mean + tmp2 * (0:(h - 1) + (1 - (1 - alpha)^n) / alpha)
  
  # Reseasonalize
  if (seasonal) {
    fcast$mean <- fcast$mean *
      rep(tail(decomp$seasonal, m), trunc(1 + h / m))[1:h]
    fcast$fitted <- fcast$fitted * decomp$seasonal
  }
  fcast$residuals <- origx - fcast$fitted
  
  # Find prediction intervals
  fcast.se <- sqrt(fcast$model$sigma2) * sqrt((0:(h - 1)) * alpha^2 + 1)
  nconf <- length(level)
  fcast$lower <- fcast$upper <- ts(matrix(NA, nrow = h, ncol = nconf))
  tsp(fcast$lower) <- tsp(fcast$upper) <- tsp(fcast$mean)
  for (i in 1:nconf) {
    zt <- -qnorm(0.5 - level[i] / 200)
    fcast$lower[, i] <- fcast$mean - zt * fcast.se
    fcast$upper[, i] <- fcast$mean + zt * fcast.se
  }
  
  # Return results
  fcast$x <- origx
  fcast$level <- level
  fcast$method <- paste0("ContextAwareTheta(", theta, ")")
  fcast$model <- list(
    alpha = alpha,
    theta = theta,
    drift = tmp2,
    sigma = fcast$model$sigma2
  )
  fcast$model$call <- match.call()
  fcast
}