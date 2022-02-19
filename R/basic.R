#' Basic forecasting ("mean", "median", "rw")
#'
#' Basic forecasting functions for multivariate time series
#'
#' @param y A multivariate time series of class \code{ts} or a matrix
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param method forecasting method, either "mean", "median", or random walk ("rw")
#' @param type_pi type of prediction interval currently, "bootstrap"
#' @param seed reproducibility seed for \code{type_pi == 'bootstrap'}
#' @param B Number of bootstrap replications for \code{type_pi == 'bootstrap'}
#'
#' @return
#' @export
#'
#' @examples
basicf <- function(y,
                   h = 5,
                   level = 95,
                   method = c("mean", "median", "rw"),
                   type_pi = c("bootstrap", "gaussian"),
                   seed = 1,
                   B = 100)
{
  if (!is.ts(y))
  {
    y <- ts(y)
  }

  method <- match.arg(method)
  stopifnot(!is.null(ncol(y)))
  n_series <- ncol(y)
  series_names <- colnames(y)
  type_pi <- match.arg(type_pi)
  freq_x <- frequency(y)
  start_x <- start(y)
  start_preds <- tsp(y)[2] + 1 / freq_x

  if (type_pi == "bootstrap")
  {
    if (method == "mean")
    {
      point_forecast <- colMeans(y)
      fcast <- t(replicate(h, point_forecast))
      resids <- NA
    }

    if (method == "median")
    {
      point_forecast <- apply(y, 2, median)
      fcast <- t(replicate(h, point_forecast))
      resids <- NA
    }

    if (method == "rw")
    {
      point_forecast <- y[nrow(y), ]
      fcast <- t(replicate(h, point_forecast))
      resids <- NA
    }
  }

  if (type_pi == "gaussian")
  {
    stop("Not implemented yet")
  }
}
