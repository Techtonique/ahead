#' Combined ets-arima-theta forecasts
#'
#' Combined ets, arima, and theta (eat) forecasting (uses \code{forecast::ets},
#' \code{forecast::auto.arima}, \code{forecast::thetaf})
#'
#' @export
#' @param y a univariate time series
#' @param h number of periods for forecasting
#' @param level confidence level for prediction intervals
#' @param method forecasting method: "E" for \code{forecast::ets};
#' "A"for \code{forecast::auto.arima}; "T" for \code{forecast::thetaf};
#'  or "EAT" for the combination of the three (default, with \code{weights})
#' @param weights weights for each method, in method \code{EAT}. Must add up to 1.
#' @param type_pi type of prediction interval:  currently ETS: "E", Auto.Arima: "A" or Theta: "T"
#' @param ... additional parameters to be passed to \code{forecast::ets},
#' \code{forecast::auto.arima}, \code{forecast::thetaf} and
#' \code{forecast::forecast}
#'
#' @details ensemble forecasts obtained from \code{forecast::ets},
#' \code{forecast::auto.arima} and \code{forecast::theta} (with weights)
#'
#' @return An object of class "forecast"; a list containing the following elements:
#'
#' \item{model}{A list containing information about the fitted model}
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{x}{The original time series}
#' \item{residuals}{Residuals from the fitted model}
#'
#' @author T. Moudiki
#'
#' @export
#'
#' @references
#'
#' Hyndman R, Athanasopoulos G, Bergmeir C, Caceres G, Chhay L,
#' O'Hara-Wild M, Petropoulos F, Razbash S, Wang E, Yasmeen F (2021).
#' forecast: Forecasting functions for time series and linear models. R
#' package version 8.14, <URL: https://pkg.robjhyndman.com/forecast/>. \cr
#'
#' Hyndman RJ, Khandakar Y (2008). 'Automatic time series forecasting: the
#' forecast package for R.' Journal of Statistical Software, 26 (3),
#' 1-22. <URL: https://www.jstatsoft.org/article/view/v027i03>.
#'
#' Assimakopoulos, V. and Nikolopoulos, K. (2000). The theta model: a
#' decomposition approach to forecasting. International Journal of
#' Forecasting 16, 521-530.
#'
#' Hyndman, R.J., and Billah, B. (2003) Unmasking the Theta method.
#' International J. Forecasting, 19, 287-290.
#'
#'
#' @examples
#'
#'require(forecast)
#'
#'\dontrun{
#'
#'print(ahead::eatf(WWWusage, method = "EAT",
#'weights = c(0.5, 0, 0.5)))
#'
#'print(ahead::eatf(WWWusage, method = "EAT"))
#'
#'
#'obj <- ahead::eatf(WWWusage, method = "EAT",
#'weights = c(0, 0.5, 0.5), h=10,
#'type_pi = "T")
#'plot(obj)
#'
#'
#'obj <- ahead::eatf(WWWusage, method = "EAT",
#'weights = c(0, 0.5, 0.5), h=10, type_pi="A")
#'plot(obj)
#'}
#'
#'
#' par(mfrow=c(3, 2))
#' plot(ahead::eatf(USAccDeaths, h=10, level=95))
#' plot(ahead::eatf(AirPassengers, h=10, level=95, type_pi = "T"))
#' plot(ahead::eatf(lynx, h=10, level=95, type_pi = "A"))
#' plot(ahead::eatf(WWWusage, h=10, level=95, type_pi = "E"))
#' plot(ahead::eatf(Nile, h=10, level=95))
#' plot(ahead::eatf(fdeaths, h=10, level=95))
#'
#'
eatf <- function(y, h = 5,
                 level = 95,
                 method = c("EAT", "E", "A", "T"),
                 weights = rep(1/3, 3),
                 type_pi = c("gaussian", "E", "A", "T"),
                 ...) {

  stopifnot(length(level) == 1)

  method <- match.arg(method)

  type_pi <- match.arg(type_pi)

  # stopifnot(length(level) == 1)
  n_y <- length(y)

  # start and frequency for returned result
  tspx <- tsp(as.ts(y))
  start_preds <- tspx[2] + 1 / tspx[3]
  freq_y <- tspx[3]

  stopifnot(sum(weights) == 1)

  if (method %in% c("E", "A", "T")) # ok, nothing to check
  {

    if (method == "E")
    {
      return(forecast::forecast(forecast::ets(y = y, ...),
                                h = h, level = level,...))
    }

    if (method == "A")
    {
     return(forecast::forecast(forecast::auto.arima(y = y, ...),
                               h = h, level = level, ...))
    }

    if (method == "T")
    {
      return(forecast::thetaf(y = y, h = h,
                              level = level, ...))
    }

  } else { #if (method == "EAT")

    stopifnot(length(weights) == 3)

    if (all(weights != 0))
    {
        obj_ets <- forecast::forecast(forecast::ets(y = y, ...),
                                      h = h, level = level, ...)
        obj_arima <- forecast::forecast(forecast::auto.arima(y = y, ...),
                                        h = h, level = level, ...)
        obj_theta <- forecast::thetaf(y = y, h = h,
                                      level = level, ...)

        fcasts <- rbind(
          E = obj_ets$mean,
          A = obj_arima$mean,
          T = obj_theta$mean)

        resids <- rbind(
          E = obj_ets$residuals,
          A = obj_arima$residuals,
          T = obj_theta$residuals)

        # not ok: forecast residuals at line 116 as in dynrm.R line 463
        out <- list()

        out$x <- y

        out$model <- list(method=method, weights=weights, type_pi=type_pi)

        out$mean <- try(ts(drop(crossprod(weights, fcasts)),
                           start = start_preds,
                           frequency = freq_y), silent = TRUE)
        if (inherits(out$mean, "try-error")) {
          out$mean <- try(ts(drop(tcrossprod(weights, fcasts)),
                           start = start_preds,
                           frequency = freq_y), silent = TRUE)
        }      
        out$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(y),
                                frequency = frequency(y))
        if (inherits(out$residuals, "try-error")) {
          out$residuals <- ts(drop(tcrossprod(weights, resids)),
                                start = start(y),
                                frequency = frequency(y))
        }                        
        out$method <- paste0("EAT(", type_pi, ")")

    } else { # if (any(weights == 0))

      if (weights[1] == 0){ # model AT

        obj_arima <- forecast::forecast(forecast::auto.arima(y = y, ...),
                                        h = h, level = level, ...)
        obj_theta <- forecast::thetaf(y = y, h = h,
                                      level = level, ...)

        fcasts <- rbind(
          E = rep(0, h),
          A = obj_arima$mean,
          T = obj_theta$mean)

        resids <- rbind(
          E = rep(0, n_y),
          A = obj_arima$residuals,
          T = obj_theta$residuals)

        out <- list()

        out$x <- y

        out$model <- list(method=method, weights=weights, type_pi=type_pi)

        out$mean <- ts(drop(crossprod(weights, fcasts)),
                           start = start_preds,
                           frequency = freq_y)
        out$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(y),
                                frequency = frequency(y))
        out$method <- paste0("AT(", type_pi, ")")
      }

      if (weights[2] == 0){ # model ET

        obj_ets <- forecast::forecast(forecast::ets(y = y, ...),
                                      h = h, level = level, ...)
        obj_theta <- forecast::thetaf(y = y, h = h,
                                      level = level, ...)

        fcasts <- rbind(
          E = obj_ets$mean,
          A = rep(0, h),
          T = obj_theta$mean)

        resids <- rbind(
          E = obj_ets$residuals,
          A = rep(0, n_y),
          T = obj_theta$residuals)


        # not ok: forecast residuals at line 116 as in dynrm.R line 463
        out <- list()

        out$x <- y

        out$model <- list(method=method, weights=weights, type_pi=type_pi)

        out$mean <- ts(drop(crossprod(weights, fcasts)),
                           start = start_preds,
                           frequency = freq_y)
        out$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(y),
                                frequency = frequency(y))
        out$method <- paste0("ET(", type_pi, ")")
      }

      if (weights[3] == 0){ # model EA

        obj_ets <- forecast::forecast(forecast::ets(y = y, ...),
                                      h = h, level = level, ...)
        obj_arima <- forecast::forecast(forecast::auto.arima(y = y, ...),
                                        h = h, level = level, ...)

        fcasts <- rbind(
          E = obj_ets$mean,
          A = obj_arima$mean,
          T = rep(0, h))

        resids <- rbind(
          E = obj_ets$residuals,
          A = obj_arima$residuals,
          T = rep(0, n_y))

        # not ok: forecast residuals at line 116 as in dynrm.R line 463
        out <- list()

        out$x <- y

        out$model <- list(method=method, weights=weights, type_pi=type_pi)

        out$mean <- ts(drop(crossprod(weights, fcasts)),
                           start = start_preds,
                           frequency = freq_y)
        out$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(y),
                                frequency = frequency(y))
        out$method <- paste0("EA(", type_pi, ")")
      }
    }

    # Compute prediction intervals
     if (type_pi == "E")
      {
        resid_fcast <- forecast::forecast(forecast::ets(out$residuals),
                                                     h = h, level=level)
      }

      if (type_pi == "A")
      {
        resid_fcast <- forecast::forecast(forecast::auto.arima(out$residuals),
                                                     h = h, level=level)
      }

      if (type_pi == "T")
      {
        resid_fcast <- forecast::thetaf(out$residuals, h = h, level=level)
      }

      if (type_pi == "gaussian")
      {
        qt_sd <- -qnorm(0.5 - level/200)*sd(out$residuals)
        rep_qt_sd <- ts(rep(qt_sd, h), start = start_preds, frequency = freq_y)

        resid_fcast <- list()
        resid_fcast$mean <- 0
        resid_fcast$lower <- -rep_qt_sd
        resid_fcast$upper <- rep_qt_sd
      }

      out$mean <-  drop(out$mean + resid_fcast$mean)
      #print("out$mean")
      #print(out$mean)
      #print("\n")
      out$lower <- drop(out$mean + resid_fcast$lower)
      #print("out$lower")
      #print(out$lower)
      #print("\n")
      out$upper <- drop(out$mean + resid_fcast$upper)

    out$level <- level

    return(structure(out, class = "forecast"))

  }

}
eatf <- compiler::cmpfun(eatf)
