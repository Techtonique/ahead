#' Plot multivariate time series forecast
#'
#' @param x result from varf or ridge2f (multivariate time series forecast)
#' @param selected_series name of the time series selected for plotting
#' @param type_graph "pi": basic prediction intervals; "dist": a distribution; "sims": the simulations
#' @param ... additional parameters to be passed to \code{plot} or \code{matplot}
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
#' obj <- ahead::ridge2f(fpp::insurance, h = 10, type_pi = "blockbootstrap", block_length=5, B = 10)
#' par(mfrow=c(1, 2))
#' plot(obj, selected_series = "Quotes", type_graph = "sims", main = "Predictive simulation for Quotes")
#' plot(obj, selected_series = "TV.advert", type_graph = "sims", main = "Predictive simulation for TV.advert")
#'
plot.mtsforecast <- function(x, selected_series,
                             type_graph = c("pi", "dist", "sims"),
                             ...)
{
  type_graph <- match.arg(type_graph)

  if (inherits(x, 'mtsforecast'))
  {
    if (type_graph == "pi")
    {
      y <- x$x[, selected_series]
      mean_fcast <- x$mean[, selected_series]
      upper_fcast <- x$upper[, selected_series]
      lower_fcast <- x$lower[, selected_series]

      start_y <- start(y)
      frequency_y <- frequency(y)

      y_mean <- ts(c(y, mean_fcast), start = start_y,
                   frequency = frequency_y)
      y_upper <- ts(c(y, upper_fcast), start = start_y,
                    frequency = frequency_y)
      y_lower <- ts(c(y, lower_fcast), start = start_y,
                    frequency = frequency_y)

      plot(y_mean, type='l',
           main=paste0("Forecasts for ", selected_series, " (", x$method, ")"),
           ylab="", ylim = c(min(c(y_upper, y_lower)),
                             max(c(y_upper, y_lower))), ...)
      lines(y_upper, col="gray60")
      lines(y_lower, col="gray60")
      lines(y_mean)
    }

    if (type_graph == "dist")
    {
      stop("Not implemented")
    }

    if (type_graph == "sims")
    {
      B <- length(x$sims)

      actual_selected_series <- x$x[, selected_series]
      series_reps <- replicate(n = B, expr = actual_selected_series)
      series_reps <- ts(series_reps,
                        start = start(actual_selected_series),
                        frequency = frequency(actual_selected_series))
      start_x <- start(x$x)
      frequency_x <- frequency(x$x)
      start_preds <- start(x$mean)

      (preds_simulations <-
          ts(
            sapply(1:B, function (i)
              x$sims[[i]][, selected_series]),
            start = start_preds,
            frequency = frequency_x
          ))

      out <- ts(
        data = rbind(series_reps, preds_simulations),
        start = start_x,
        frequency = frequency_x
      )

      time_inputs <- as.numeric(time(x$x))
      input_series <- x$x[, selected_series]

      matplot(
        x = as.numeric(time(out)),
        y = out,
        xlab = 'time',
        ylab = selected_series,
        type = 'l',
        lwd = 2,
        ...
      )
      lines(x = time_inputs, y = input_series, lwd = 2)
    }


  } else {
   stop("Method not implemented for this type of object")
  }
}
