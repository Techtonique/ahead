
#' Loess forecasting
#'
#' @param y
#' @param h
#' @param level
#' @param span
#' @param degree
#' @param b
#' @param B
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#'
#' plot(loessf(Nile, h=20, level=95))
#'
loessf <- function(y, h = 5, level = 95,
                   span = 0.75, degree = 2,
                   b = NULL, B = 250,
                   seed = 123,
                   type_boot = c("pkgforecast", "pkgahead"))
{
  # adapted from forecast:::bld.mbb.bootstrap ---
  freq_y <- frequency(y)
  if (length(y) <= 2 * freq_y)
    freq_y <- 1L

  if (is.null(b)) {
    b <- ifelse(freq_y > 1, 2 * freq_y, min(8, floor(length(y)/2)))
  }

  type_boot <- match.arg(type_boot)

  # start and frequency for returned result
  tspy <- tsp(as.ts(y))
  start_preds <- tspy[2] + 1 / tspy[3]

  # Adjust
  input_times <- time(y)
  fit_loess <- stats::loess(value ~ time,
               data = cbind.data.frame(time = input_times,
                                       value = y),
               span = span,
               degree = degree,
               control = loess.control(surface = "direct"))
  resids <- ts(stats::residuals(fit_loess), start = start(y),
               frequency = freq_y)
  fitted_values <- ts(stats::fitted(fit_loess), start = start(y),
                      frequency = freq_y)

  # point forecast

  time_y <- time(y)

  diff_time_y <- diff(time_y)[1]

  time_oos <- seq(time_y[length(time_y)] + diff_time_y,
                  by = diff_time_y,
                  length.out = h)

  fcast <- ts(as.numeric(predict(fit_loess,
                                 data.frame(time = time_oos))),
              start = start_preds, frequency = freq_y)

  simulations <- ts(matrix(NA, nrow = h, ncol = B),
             start = start_preds,
             frequency = freq_y)

  if (type_boot == "pkgahead")
  {
    for (i in 1:B)
    {
      # sampling from the residuals
      bootstrapped_residuals <- ts(drop(mbb(r = matrix(resids, ncol = 1),
                                            n = h, b = b, seed=100*i+3)),
                                   start = start_preds,
                                   frequency = freq_y)

      simulations[, i] <- ts(fcast + bootstrapped_residuals,
                             start = start_preds,
                             frequency = freq_y)
    }
  }

  if (type_boot == "pkgforecast")
  {
    for (i in 1:B)
    {
    # sampling from the residuals
    bootstrapped_residuals <- ts(mbb2(x = resids,
                                      window_size = b),
                                 start = start_preds,
                                 frequency = freq_y)

    simulations[, i] <- ts(fcast + bootstrapped_residuals,
                           start = start_preds,
                           frequency = freq_y)
    }
  }

  preds_upper <-  preds_lower <- rep(0, h)

  preds_mean <- ts(apply(simulations, 1, median),
                   start = start_preds,
                   frequency = freq_y)

  preds_upper <-
    ts(apply(simulations, 1, function(x)
      quantile(x, probs = 1 - (1 - level / 100) / 2)),
      start = start_preds,
      frequency = freq_y)

  preds_lower <-
    ts(apply(simulations, 1, function(x)
      quantile(x, probs = (1 - level / 100) / 2)),
      start = start_preds,
      frequency = freq_y)

  out <- list(mean = preds_mean,
              lower = preds_lower,
              upper = preds_upper,
              level = level,
              fitted = fitted_values,
              residuals = resids,
              model = fit_loess,
              method = "loessf",
              sims = simulations,
              x = y)

  # 'forecast' must be in the environment
  return(structure(out, class = "forecast"))
}
