#' Plot multivariate time series forecast
#'
#' @param x result from varf or ridge2f (multivariate time series forecast)
#' @param selected_series name of the time series selected for plotting
#' @param ... additional parameters to be passed to \code{plot}
#'
#' @export
#'
#' @examples
#'
#' require(fpp)
#'
#' fit_obj_VAR <- ahead::varf(fpp::insurance, lags = 2,
#' h = 10, level = 95)
#'
#' fit_obj_ridge2 <- ahead::ridge2f(fpp::insurance, lags = 2,
#' h = 10, level = 95)
#'
#' par(mfrow=c(2, 2))
#' plot(fit_obj_VAR, "Quotes")
#' plot(fit_obj_VAR, "TV.advert")
#' plot(fit_obj_ridge2, "Quotes")
#' plot(fit_obj_ridge2, "TV.advert")
#'
plot.mtsforecast <- function(x, selected_series, ...)
{
  if (methods::is(x, 'mtsforecast'))
  {
    y <- x$x[, selected_series]
    mean_fcast <- x$mean[, selected_series]
    upper_fcast <- x$upper[, selected_series]
    lower_fcast <- x$lower[, selected_series]

    y_mean <- ts(c(y, mean_fcast), start = start(y),
                 frequency = frequency(y))
    y_upper <- ts(c(y, upper_fcast), start = start(y),
                  frequency = frequency(y))
    y_lower <- ts(c(y, lower_fcast), start = start(y),
                  frequency = frequency(y))

    plot(y_mean, type='l',
         main=paste0("Forecasts for ", selected_series, " (", x$method, ")"),
         ylab="", ylim = c(min(c(y_upper, y_lower)),
         max(c(y_upper, y_lower))), ...)
    lines(y_upper, col="gray60")
    lines(y_lower, col="gray60")
    lines(y_mean)
  } else {
   stop("Method not implemented for this type of object")
  }
}
