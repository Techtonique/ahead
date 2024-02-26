#' ARMA(1, 1)-GARCH(1, 1) forecasting (with simulation)
#'
#' @param y a univariate time series
#' @param h number of periods for forecasting
#' @param level confidence level for prediction intervals
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
#' log_returns <- ts(log(y[-1]/y[-length(y)]))
#'
#' # require(forecast)
#' # z <- ahead::armagarchf(y=log_returns, h=200)
#' # plot(z)
#'
#'
armagarchf <- function(y,
                       h = 5,
                       level = 95,
                       B = 250,
                       cl = 1L,
                       dist = c("student", "gaussian"),
                       seed = 123) {

  if(is_package_available("forecast") == FALSE)
    utils::install.packages("forecast",
                            repos = c(CRAN = "https://cloud.r-project.org"),
                            dependencies = TRUE)

  if(is_package_available("fBasics") == FALSE)
    utils::install.packages("fBasics",
                            repos = c(CRAN = "https://cloud.r-project.org"),
                            dependencies = TRUE)

  if(is_package_available("fGarch") == FALSE)
    utils::install.packages("fGarch",
                            repos = c(CRAN = "https://cloud.r-project.org"),
                            dependencies = TRUE)

  dist <- match.arg(dist)

  # Fit ARIMA to original series
  obj_arma <- stats::arima(x = y, order = c(1L, 0L, 1L))

  # In sample residuals
  eps <- stats::residuals(obj_arma)

  eps_prev <- eps[length(eps)]

  # Fit GARCH to residuals
  base::suppressWarnings(
    obj_garch <- fGarch::garchFit(
      formula =  ~ garch(1, 1),
      data = eps,
      include.mean = FALSE,
      trace = FALSE
    )
  )
  omega <- obj_garch@fit$coef["omega"]
  alpha <- obj_garch@fit$coef["alpha1"]
  beta <- obj_garch@fit$coef["beta1"]


  # Forecasting using simulation
  if (cl <= 1L)
  {
    res <-
      sapply(1:B, function(x)
        stats::arima.sim(
          model = list(ar = obj_arma$coef["ar1"],
                       ma = obj_arma$coef["ma1"]),
          n = h,
          start.innov = eps,
          innov = forecast_innovs(
            eps_prev = eps_prev,
            h =
              h,
            omega = omega,
            alpha =
              alpha,
            beta = beta,
            seed =
              x + 10 * seed,
            dist = dist
          )
        ))
  } else {
    stopifnot(is.numeric(cl))
    cluster <- parallel::makeCluster(getOption("cl.cores", cl))
    res <-
      parallel::parSapply(
        cl = cluster,
        X = 1:B,
        FUN = function(x)
          stats::arima.sim(
            model = list(ar = obj_arma$coef["ar1"],
                         ma = obj_arma$coef["ma1"]),
            n = h,
            start.innov = eps,
            innov = forecast_innovs(
              eps_prev = eps_prev,
              h =
                h,
              omega = omega,
              alpha =
                alpha,
              beta = beta,
              seed =
                x + 10 * seed,
              dist = dist
            )
          )
      )
    parallel::stopCluster(cluster)
  }


  # return
  tspx <- tsp(y)
  start_preds <- tspx[2] + 1 / tspx[3]
  freq_y <- tspx[3]
  ans <- list()
  ans$x <- y

  ans$mean <- ts(apply(res, 1, median),
                 start = start_preds,
                 frequency = freq_y)

  ans$level <- level
  level_100 <- level / 100
  qts <-
    apply(res, 1, function(x)
      quantile(x, probs = c((1 - level_100) / 2,
                            1 - (1 - level_100) /
                              2)))
  ans$lower <- ts(qts[1,], start = start_preds, frequency = freq_y)
  ans$upper <- ts(qts[2,], start = start_preds, frequency = freq_y)

  ans$method <- "ARMA(1, 1) - GARCH(1, 1)"
  ans$model <- list(
    h = h,
    level = level,
    B = B,
    seed = seed
  )
  ans$sims <- res

  return(structure(ans, class = "forecast"))
}
armagarchf <- compiler::cmpfun(armagarchf)


# Utils function for armagarchf
forecast_innovs <- function(eps_prev,
                            omega,
                            alpha,
                            beta,
                            df = 3,
                            h = 5,
                            seed = 123,
                            dist = c("student", "gaussian"))
{
  sigma_prev <- omega / (1 - alpha - beta)

  eps <- rep(NA, h)

  set.seed(seed)

  dist <- match.arg(dist)

  if (dist == "student")
  {
    rts <- rt(n = h, df = df)

    eps <- forecast_innovs_loop_cpp(eps, rts,
                                    eps_prev,
                                    omega, alpha, beta,
                                    df, h)
  } else {
    rn <- rnorm(n = h)

    eps <- forecast_innovs_loop_cpp2(eps, rn,
                                     eps_prev,
                                     omega, alpha, beta,
                                     df, h)
  }



  return(drop(eps))
}
