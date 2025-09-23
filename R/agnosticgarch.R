#' ANY MODEL+GARCH(1, 1) forecasting 
#'
#' @param y a univariate time series
#' @param h number of periods for forecasting
#' @param level confidence level for prediction intervals
#' @param FUN forecasting function for the main model; 
#' default \code{ahead::dynrmf}
#' @param B number of simulations for `arima.sim`
#' @param cl an integer; the number of clusters for parallel execution
#' @param dist distribution of innovations ("student" or "gaussian")
#' @param seed reproducibility seed
#'
#' @return An object of class "forecast"; a list containing the following elements:
#'
#' \item{model}{A list containing information about the fitted model}
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{x}{The original time series}
#' \item{sims}{Simulations of ARMA(1, 1)-GARCH(1, 1)}
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
                           dist = c("student", "gaussian"),
                           seed = 123, 
                           ...) {
  
  dist <- match.arg(dist)
  # Fit ARIMA to original series
  obj_mean <- ahead::genericforecast(FUN=FUN, y=y, h=h, ...)
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
    preds_vol <- fGarch::predict(obj_garch, n.ahead=h, conf=level/100, plot=TRUE); 
    dev.off();
    ans$mean <- ts(as.numeric(preds_main + preds_vol$meanForecast), 
                   start = start_preds, frequency = freq_y)
    ans$lower <- ts(as.numeric(preds_main + preds_vol$lowerInterval), 
                    start = start_preds, frequency = freq_y)
    ans$upper <- ts(as.numeric(preds_main + preds_vol$upperInterval), 
                    start = start_preds, frequency = freq_y)
    ans$method <- paste0(obj_mean$method, "+GARCH(1, 1)")
  return(structure(ans, class = "forecast"))
}