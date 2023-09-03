

#' Loess forecasting
#'
#' @param y A numeric vector or time series of class \code{ts}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param span the parameter which controls the degree of smoothing
#' @param degree the degree of the polynomials to be used, normally 1 or 2. (Degree 0 is also allowed, but see \code{stats::loess})
#' @param type_pi Type of prediction interval currently (independent) "bootstrap", (circular) "blockbootstrap", or "movingblockbootstrap"
#' @param b block length for circular block bootstrap
#' @param B number of bootstrap replications
#' @param type_aggregation Type of aggregation, ONLY for bootstrapping; either "mean" or "median"
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
#' \item{residuals}{Residuals from the fitted model}
#' \item{sims}{Model simulations for bootstrapping}
#'
#' @export
#'
#' @examples
#'
#' par(mfrow = c(3, 1))
#'
#' plot(loessf(Nile, h=20, level=95, B=10))
#'
#' plot(loessf(Nile, h=20, level=95, B=10,
#'      type_pi = "blockbootstrap"))
#'
#' plot(loessf(Nile, h=20, level=95, B=10,
#'      type_pi = "movingblockbootstrap"))
#'
loessf <- function(y,
                   h = 5,
                   level = 95,
                   span = 0.75,
                   degree = 2,
                   type_pi = c("bootstrap",
                               "blockbootstrap",
                               "movingblockbootstrap"),
                   b = NULL,
                   B = 250,
                   type_aggregation = c("mean", "median"),
                   seed = 123)
{
  freq_y <- frequency(y)
  if (length(y) <= 2 * freq_y)
    freq_y <- 1L

  if (is.null(b)) {
    b <- ifelse(freq_y > 1, 2 * freq_y, min(8, floor(length(y) / 2)))
  }

  type_pi <- match.arg(type_pi)
  type_aggregation <- match.arg(type_aggregation)

  # start and frequency for returned result
  tspy <- tsp(as.ts(y))
  start_preds <- tspy[2] + 1 / tspy[3]

  # Adjust
  input_times <- stats::time(y)
  fit_loess <- stats::loess(
    value ~ time,
    data = cbind.data.frame(time = input_times,
                            value = y),
    span = span,
    degree = degree,
    control = stats::loess.control(surface = "direct")
  )
  resids <- ts(stats::residuals(fit_loess),
               start = start(y),
               frequency = freq_y)
  fitted_values <- ts(stats::fitted(fit_loess),
                      start = start(y),
                      frequency = freq_y)

  # point forecast

  time_y <- time(y)

  diff_time_y <- diff(time_y)[1]

  time_oos <- seq(time_y[length(time_y)] + diff_time_y,
                  by = diff_time_y,
                  length.out = h)

  fcast <- ts(as.numeric(predict(fit_loess,
                                 data.frame(time = time_oos))),
              start = start_preds,
              frequency = freq_y)

  simulations <- ts(matrix(NA, nrow = h, ncol = B),
                    start = start_preds,
                    frequency = freq_y)

  for (i in 1:B)
  {
    set.seed(100 * i + 3)
    # sampling from the residuals
    idx <- sample.int(n = length(resids),
                      size = h,
                      replace = TRUE)

    bootstrapped_residuals <- ts(switch(
      type_pi,
      bootstrap = matrix(resids, ncol = 1L)[idx,],
      blockbootstrap = drop(
        mbb(
          r = matrix(resids, ncol = 1L),
          n = h,
          b = b,
          seed = 100 * i + 3,
          return_indices =
            FALSE
        )),
        movingblockbootstrap = drop(
          mbb2(
            r = matrix(resids, ncol = 1L),
            n = h,
            b = b,
            seed = 100 * i + 3,
            return_indices =
              FALSE
          ))
      ),
    start = start_preds,
    frequency = freq_y)

    simulations[, i] <- ts(fcast + bootstrapped_residuals,
                           start = start_preds,
                           frequency = freq_y)
  }

  preds_upper <-  preds_lower <- rep(0, h)

  preds_mean <- ts(switch(type_aggregation,
                          mean = rowMeans(simulations),
                          median = apply(simulations, 1, median)),
                   start = start_preds,
                   frequency = freq_y)

  preds_upper <-
    ts(
      apply(simulations, 1, function(x)
        quantile(x, probs = 1 - (1 - level / 100) / 2)),
      start = start_preds,
      frequency = freq_y
    )

  preds_lower <-
    ts(
      apply(simulations, 1, function(x)
        quantile(x, probs = (1 - level / 100) / 2)),
      start = start_preds,
      frequency = freq_y
    )

  out <- list(
    mean = preds_mean,
    lower = preds_lower,
    upper = preds_upper,
    level = level,
    fitted = fitted_values,
    residuals = resids,
    model = fit_loess,
    method = "loessf",
    sims = simulations,
    x = y
  )

  # 'forecast' must be in the environment
  return(structure(out, class = "forecast"))
}
