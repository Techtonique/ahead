#' ANY MODEL+GARCH(1, 1) forecasting 
#'
#' @param y a univariate time series
#' @param h number of periods for forecasting
#' @param level confidence level for prediction intervals
#' @param FUN forecasting function for the main model; 
#' default \code{ahead::dynrmf}
#' @param seed reproducibility seed
#' @param ... 
#' 
#' @return An object of class "forecast"; a list containing the following elements:
#'
#' \item{model}{A list containing information about the fitted model}
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{x}{The original time series}
#' \item{sims}{Simulations of ANYMODEL+GARCH(1, 1)}
#'
#' @author T. Moudiki
#'
#' @export
#'
#' @examples
#'
#' y <- datasets::EuStockMarkets[ , "DAX"]
#'
#' # require(forecast)
#' # z <- ahead::agnosticgarchf(y=y, h=20)
#' # plot(z)
#'
agnosticgarchf <- function(y,
                           h = 5,
                           level = 95,
                           FUN=forecast::auto.arima,
                           seed = 123, 
                           ...) {
  
  # Fit ARIMA to original series
  obj_mean <- try(ahead::genericforecast(FUN=FUN, y=y, h=h, level=95, ...), 
                  silent=TRUE)
  if (inherits(obj_mean, "try-error"))
  {
    obj_mean <- FUN(y=y, h=h, level=95, ...)
  }
  # In sample residuals
  eps <- residuals(obj_mean)
  eps_prev <- eps[length(eps)]
  # Fit GARCH to residuals
  base::suppressWarnings(
    obj_garch <- fGarch::garchFit(
      formula =  ~ garch(1, 1),
      data = eps,
      include.mean = FALSE,
      trace = FALSE, 
      ...
    )
  )
  
    tspx <- tsp(y)
    start_preds <- tspx[2] + 1 / tspx[3]
    freq_y <- tspx[3]
    ans <- list()
    ans$x <- y
    ans$level <- level
    ans$model <- list(obj_mean=obj_mean, 
                      obj_garch=obj_garch)
    preds_main <- obj_mean$mean
    # Get distribution info
    dist_ <- obj_garch@fit$params$cond.dist
    alpha <- 1 - level/100  # for 95% PI
    if (dist_ %in% c("norm", "QMLE")) {
      z <- qnorm(1 - alpha/2)
    } else if (dist_ == "std") {
      nu <- coef(fit)["nu"]
      if (is.na(nu)) nu <- fit@fit$par["nu"]  # safeguard
      z <- qt(1 - alpha/2, df = nu)
    } else if (dist_ == "ged") {
      # GED is trickier; fGarch doesn't export qged easily
      # Approximate with normal or use external package like `fBasics`
      z <- qnorm(1 - alpha/2)  # common practical choice
    } else {
      warning("Unknown distribution; using normal approximation.")
      z <- qnorm(1 - alpha/2)
    }
    preds_vol <- fGarch::predict(obj_garch, n.ahead = h, conf = level/100, plot = FALSE)
    # Compute intervals
    #fitted_values <- try(fitted(obj_mean), silent=TRUE) + try(fitted(obj_garch), silent=TRUE)
    #resids <- try(residuals(obj_mean), silent=TRUE) + try(residuals(obj_garch), silent=TRUE)
    #ans$fitted <- fitted_values
    #ans$residuals <- resids
    se <- preds_vol$meanError
    lower <- preds_vol$meanForecast - z * se
    upper <- preds_vol$meanForecast + z * se
    ans$mean <- ts(as.numeric(preds_main + preds_vol$meanForecast), 
                   start = start_preds, frequency = freq_y)
    ans$lower <- ts(as.numeric(preds_main + lower), 
                    start = start_preds, frequency = freq_y)
    ans$upper <- ts(as.numeric(preds_main + upper), 
                    start = start_preds, frequency = freq_y)
    ans$method <- paste0(obj_mean$method, "+GARCH(1, 1)")
  return(structure(ans, class = "forecast"))
}