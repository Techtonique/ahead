#' Compute First-Order Sensitivity Effects for Dynamic Regression Forecasts
#'
#' This function computes the sensitivity of a forecasted time series (`y`) 
#' to small perturbations in an external regressor (`xreg`) using the `dynrmf` function.
#' 
#' @param y A univariate time series object to forecast.
#' @param xreg A time series of external regressors (same length as `y`).
#' @param h Forecast horizon (length of test set). If NULL, uses length of xreg.
#' @param level Confidence level for `dynrmf` (default 99).
#' @param zero Small positive constant to avoid division by zero (default 1e-4).
#' @param ... Additional parameters to be passed to \code{\link{ahead::dynrmf}}
#' 
#' @return A list containing:
#' \describe{
#'   \item{effects_mean}{Sensitivity effects on the mean forecast.}
#'   \item{forecast_central}{Forecast without xreg.}
#'   \item{forecast_with_xreg}{Forecast with xreg.}
#'   \item{forecast_plus}{Forecast with positive perturbation on xreg.}
#'   \item{forecast_minus}{Forecast with negative perturbation on xreg.}
#'   \item{h}{Forecast horizon used.}
#'   \item{xreg_test}{Test set of xreg used in forecasting.}
#'   \item{hh}{Perturbation applied to xreg.}
#'   \item{split}{List with training and testing sets of y.}
#' }
#' 
#' @export
dynrmf_sensi <- function(y, xreg, h = NULL, level = 99, zero = 1e-4, ...) {
  
  n <- length(y)
  if (is.null(h)) {
    h <- n
  }
  if (h > n) {
    warning("Forecast horizon h exceeds series length. Using full series.")
    h <- n
  }
  
  # Define training and test sets
  y_train <- window(y, end = time(y)[n - h])
  y_test  <- window(y, start = time(y)[n - h + 1])
  
  xreg_train <- window(xreg, end = time(xreg)[n - h])
  xreg_test  <- window(xreg, start = time(xreg)[n - h + 1])
  
  eps_factor <- zero^(1/3)
  
  # Base forecast without xreg
  central_forecast <- ahead::dynrmf(y = y_train,
                                    level = level,
                                    h = h)
  
  # Forecast with actual xreg
  central_forecast_xreg <- ahead::dynrmf(y = y_train,
                                         xreg_fit = xreg_train,
                                         xreg_predict = xreg_test,
                                         level = level,
                                         h = h, 
                                         ...)
  
  # Create perturbations
  cond <- abs(xreg_test) > zero
  hh <- eps_factor * xreg_test * cond + zero * (!cond)
  
  # Forecast with perturbed xreg
  forecast_plus_hh <- ahead::dynrmf(y = y_train,
                                    xreg_fit = xreg_train,
                                    xreg_predict = xreg_test + hh,
                                    level = level,
                                    h = h)
  
  forecast_minus_hh <- ahead::dynrmf(y = y_train,
                                     xreg_fit = xreg_train,
                                     xreg_predict = xreg_test - hh,
                                     level = level,
                                     h = h)
  
  # Compute sensitivity effects
  effects_mean <- (forecast_plus_hh$mean - forecast_minus_hh$mean) / (2 * hh)
  
  # Return results
  return(list(
    effects_mean = effects_mean,
    forecast_central = central_forecast,
    forecast_with_xreg = central_forecast_xreg,
    forecast_plus = forecast_plus_hh,
    forecast_minus = forecast_minus_hh,
    h = h,
    xreg_test = xreg_test,
    hh = hh,
    split = list(training = y_train, testing = y_test)
  ))
}

#' Plot First-Order Sensitivity Effects
#'
#' This function plots the computed sensitivity effects over the forecast horizon.
#' 
#' @param sensitivity_results Output list from `dynrmf_sensi`.
#' @param title Character. Plot title.
#' @param y_label Character. Label for y-axis.
#' @param x_label Character. Label for x-axis (default "Forecast Horizon").
#' 
#' @return ggplot2 object displaying the sensitivity effects.
#' 
#' @examples
#' \dontrun{
#' plot1 <- plot_dynrmf_sensitivity(res, 
#'                           title = "Sensitivity of Consumption to Income",
#'                           y_label = "Effect (ΔConsumption / ΔIncome)")
#' print(plot1)
#' }
#' 
#' @export
plot_dynrmf_sensitivity <- function(sensitivity_results, 
                             title = "First-Order Sensitivity Effect",
                             y_label = "Effect (ΔY / ΔX)",
                             x_label = "Forecast Horizon") {
  
  h <- sensitivity_results$h
  horizon <- 1:h
  effects_mean <- as.numeric(sensitivity_results$effects_mean)
  
  if (length(effects_mean) != h) {
    warning(paste("Length of effects_mean (", length(effects_mean), 
                  ") doesn't match h (", h, "). Truncating."))
    effects_mean <- effects_mean[1:min(length(effects_mean), h)]
    horizon <- 1:length(effects_mean)
    h <- length(effects_mean)
  }
  
  df_main <- data.frame(
    Horizon = horizon,
    Effect = effects_mean
  )
  
  ggplot(df_main, aes(x = Horizon, y = Effect)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 14)
}

#' # ============================================
#' # Example Usage
#' # ============================================
#' 
#' # Example 1: US Consumption vs Income
#' sensitivity_results_auto <- ahead::dynrmf_sensi(
#'   y = fpp2::uschange[, "Consumption"],
#'   xreg = fpp2::uschange[, "Income"],
#'   h = 10
#' )
#' 
#' plot1 <- ahead::plot_dynrmf_sensitivity(sensitivity_results_auto, 
#'                           title = "Sensitivity of Consumption to Income",
#'                           y_label = "Effect (ΔConsumption / ΔIncome)")
#' print(plot1)
#' 
#' #' Example 2: TV Advertising vs Insurance Quotes
#' sensitivity_results_tv <- ahead::dynrmf_sensi(
#'   y = fpp2::insurance[, "Quotes"],
#'   xreg = fpp2::insurance[, "TV.advert"],
#'   h = 8
#' )
#' 
#' plot3 <- plot_dynrmf_sensitivity(sensitivity_results_tv,
#'                           title = "Sensitivity of Insurance Quotes to TV Advertising",
#'                           y_label = "Effect (ΔQuotes / ΔTV.advert)")
#' print(plot3)