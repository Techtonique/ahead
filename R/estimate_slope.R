#' Estimate local slope of a fitted model's response surface with respect to time
#'
#' Fits a user-supplied model to a univariate time series (using time index
#' and an auxiliary random covariate as predictors), then numerically
#' differentiates the fitted model's predictions with respect to time using a
#' central finite-difference (perturbation) scheme. Slopes are computed for
#' every historical time point as well as for \code{h} future time points,
#' making this useful for characterizing local trend/momentum in the series
#' both in-sample and over a forecast horizon.
#'
#' @details
#' The function builds a data frame with columns \code{y} (the response,
#' padded with \code{NA} for the \code{h} forecast periods), \code{t} (an
#' integer time index running from \code{0} to \code{n + h - 1}), and
#' \code{z} (an i.i.d. standard normal covariate used only to satisfy models
#' that require more than one predictor). The model is fit on the first
#' \code{n} rows (the historical data) via \code{fit_func}, trying in order:
#' \itemize{
#'   \item the formula interface \code{fit_func(y ~ ., data = df_train, ...)};
#'   \item the matrix/vector interface \code{fit_func(x = cbind(t, z), y = y, ...)};
#'   \item the positional interface \code{fit_func(cbind(t, z), y, ...)}.
#' }
#' If none of these succeed, the function stops with an error.
#'
#' For every time point (historical and future), an adaptive step size
#' \code{h_eps = max(zero^(1/3) * |t|, zero)} (with \code{zero = 1e-4}) is
#' used to perturb \code{t} up and down while holding \code{z} fixed. The
#' fitted model is used to predict the response at \code{t + h_eps} and
#' \code{t - h_eps}, and the central-difference slope
#' \code{(fx_plus - fx_minus) / (2 * h_eps)} is returned for each time point.
#'
#' A fixed random seed (\code{123}) is set internally so that the random
#' covariate \code{z} — and therefore the resulting slopes — are reproducible
#' across calls.
#'
#' @param fit_func A function used to fit the model. It is tried successively
#'   with a formula interface (\code{y ~ ., data = ...}), an \code{x}/\code{y}
#'   interface, and a positional \code{(x, y)} interface, so it should support
#'   at least one of these calling conventions (e.g. \code{lm}, \code{randomForest},
#'   \code{glmnet}, etc.).
#' @param predict_func A function used to generate predictions from the
#'   object returned by \code{fit_func}. It is called as
#'   \code{predict_func(object, newdata)} and, if that fails, as
#'   \code{predict_func(object, as.matrix(newdata))}. Must return a numeric
#'   vector of predictions.
#' @param y A numeric vector containing the univariate time series of
#'   historical observations.
#' @param h Integer; forecasting horizon, i.e. the number of future time
#'   points (beyond the length of \code{y}) for which slopes should also be
#'   computed. Defaults to \code{0} (in-sample slopes only).
#' @param ... Additional arguments passed on to \code{fit_func}.
#'
#' @return A numeric vector of length \code{length(y) + h} containing the
#'   estimated local slope (numerical derivative of the fitted model's
#'   predictions with respect to time) at each historical and future time
#'   point, in chronological order.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' y <- cumsum(rnorm(50))
#' slopes <- estimate_theta_slope(
#'   fit_func = lm,
#'   predict_func = predict,
#'   y = y,
#'   h = 10
#' )
#' plot(slopes, type = "l")
#' }
#'
#' @export
estimate_theta_slope <- function(fit_func,
                                 predict_func,
                                 y,
                                 h = 0,
                                 ...) {
  set.seed(123)
  n <- length(y)
  
  # Extend time range to include forecast horizon
  time_idx <- 0:(n + h - 1)
  n_total <- length(time_idx)
  
  zero <- 1e-4
  eps_factor <- zero^(1 / 3)
  random_covariate <- rnorm(n_total)
  
  # Base data - only first n points have y values
  df <- data.frame(
    y = c(y, rep(NA, h)),
    t = time_idx, 
    z = random_covariate
  )
  
  # Fit model on historical data only (first n rows)
  df_train <- df[1:n, ]
  tmp2 <- try(fit_func(y ~ ., data = df_train, ...), silent = TRUE)
  if (inherits(tmp2, "try-error")) {
    tmp2 <- try(fit_func(x = cbind(df_train$t, df_train$z), 
                         y = df_train$y, ...), silent = TRUE)
    if (inherits(tmp2, "try-error")) {
      tmp2 <- try(fit_func(cbind(df_train$t, df_train$z), 
                           df_train$y, ...), silent = TRUE)
      if (inherits(tmp2, "try-error")) {
        stop("unable to fit the model")
      }
    }
  }
  
  # Adaptive step for ALL time points (historical + future)
  h_eps <- pmax(eps_factor * abs(time_idx), zero)
  double_h <- 2 * h_eps
  
  # Perturbations for all time points
  t_plus  <- time_idx + h_eps
  t_minus <- time_idx - h_eps
  df_plus  <- data.frame(t = t_plus, z = random_covariate)
  df_minus <- data.frame(t = t_minus, z = random_covariate)
  
  # Predict at perturbed points
  fx_plus <- try(as.numeric(predict_func(tmp2, df_plus)), silent = TRUE)
  if (inherits(fx_plus, "try-error"))
    fx_plus <- as.numeric(predict_func(tmp2, as.matrix(df_plus)))
  fx_minus <- try(as.numeric(predict_func(tmp2, df_minus)), silent = TRUE)
  if (inherits(fx_minus, "try-error"))
    fx_minus <- as.numeric(predict_func(tmp2, as.matrix(df_minus)))
  
  # Slopes for ALL time points (historical + future)
  slope <- (fx_plus - fx_minus) / double_h
  
  return(slope)
}

