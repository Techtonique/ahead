#' Vector Autoregressive model (adapted from vars::VAR)
#'
#' Vector Autoregressive model adapted from vars::VAR (only for benchmarking)
#'
#' @param y A multivariate time series of class \code{ts}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags
#' @param type_VAR Type of deterministic regressors to include.
#' @param ... Additional parameters to be passed to vars::VAR.
#'
#' @return An object of class "mtsforecast"; a list containing the following elements:
#'
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{x}{The original time series}
#' \item{residuals}{Residuals from the fitted model}
#'
#' @author T. Moudiki
#'
#' @references
#'
#' Bernhard Pfaff (2008). VAR, SVAR and SVEC Models: Implementation
#' Within R Package vars. Journal of Statistical Software 27(4). URL
#' http://www.jstatsoft.org/v27/i04/. \cr
#'
#' Pfaff, B. (2008) Analysis of Integrated and Cointegrated Time Series
#' with R. Second Edition. Springer, New York. ISBN 0-387-27960-1
#'
#' @export
#'
#' @examples
#'
#' require(fpp)
#'
#' print(varf(fpp::insurance, lags=2, h=10))
#'
#' res <- varf(fpp::usconsumption, h=20, lags=2)
#'
#' par(mfrow=c(1, 2))
#' plot(res, "consumption")
#' plot(res, "income")
#'
varf <- function(y, h = 5,
                 level = 95,
                 lags = 1,
                 type_VAR = c("const", "trend",
                              "both", "none"),
                 ...)
{
  if(is_package_available("vars") == FALSE)
    install.packages("vars",
                     repos = c(CRAN = "https://cloud.r-project.org"))

  if (!is.ts(y))
  {
    y <- ts(y)
  }
  freq_y <- frequency(y)
  start_y <- start(y)
  start_preds <- tsp(y)[2] + 1/freq_y

  type_VAR <- match.arg(type_VAR)

  # Fitting a VAR model to multiple time series
  fit_obj <- fit_var_mts(y, lags = lags,
                         type_VAR = type_VAR,
                         ...)

  # Forecast from fit_obj
  preds <- fcast_var_mts(fit_obj$fit_obj,
                         h = h, level = level)

  out <- list(mean = ts(sapply(preds$fcst, function(x) x[,"fcst"]),
                        start = start_preds, frequency = freq_y),
              lower = ts(sapply(preds$fcst, function(x) x[,"lower"]),
                         start = start_preds, frequency = freq_y),
              upper = ts(sapply(preds$fcst, function(x) x[,"upper"]),
                         start = start_preds, frequency = freq_y),
              x = y,
              level = level,
              method = "VAR",
              residuals = residuals(fit_obj$fit_obj))

  # cat("Point Forecast", "\n")
  # print(out$mean)
  # cat("\n")
  #
  # cat("Lo", level, "\n")
  # print(out$lower)
  # cat("\n")
  #
  # cat("Hi", level, "\n")
  # print(out$upper)
  # cat("\n")

  return(structure(out, class = "mtsforecast"))
}
varf <- compiler::cmpfun(varf)


# Fitting a VAR model to multiple time series
fit_var_mts <- function(x,
                        lags = 1,
                        type_VAR = c("const", "trend",
                                     "both", "none"),
                        ...)
  # for penalization == "none" only
{
  series_names <- colnames(x)

  # unrestricted VAR algo from package 'vars'
    type_VAR <- match.arg(type_VAR)
    fit_obj <- vars::VAR(y = x, p = lags, type = type_VAR,
                         ...)

    return(
      list(
        fit_obj = fit_obj,
        series_names = series_names,
        class_res = "VAR"
      )
    )
}


# Forecasting a VAR model
fcast_var_mts <- function(fit_obj,
                          h = 5,
                          level = 95)
{
   return(predict(fit_obj, n.ahead = h,
                     ci = level / 100))
}
