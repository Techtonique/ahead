#' Fit and forecast for benchmarking purposes
#'
#' @param y A univariate time series of class \code{ts}
#' @param h Forecasting horizon (default is \code{NULL}, in that case, \code{pct_train}
#' and \code{pct_calibration} are used)
#' @param pct_train Percentage of data in the training set, when \code{h} is \code{NULL}
#' @param pct_calibration Percentage of data in the calibration set for conformal prediction
#' @param method For now "thetaf" (default), "arima", "ets", "tbats", "tslm", "dynrmf" (from ahead),
#' "ridge2f" (from ahead), "naive", "snaive"
#' @param level Confidence level for prediction intervals in %, default is 95
#' @param B Number of bootstrap replications or number of simulations
#' (yes, 'B' is unfortunate)
#' @param seed Reproducibility seed
#' @param graph Plot or not?
#' @param conformalize Calibrate or not?
#' @param type_calibration "splitconformal" (default conformal method), "cv1" (do not use),
#' "loocv" (do not use)
#' @param gap length of training set for loocv conformal (do not use)
#' @param agg "mean" or "median" (aggregation method)
#' @param vol "constant" or "garch" (type of volatility modeling for calibrated residuals)
#' @param type_sim "kde", "surrogate", "bootstrap" (type of simulation for calibrated residuals)
#' @param ... additional parameters
#'
#' @return an object of class 'forecast' with additional information
#' @export
#'
#' @examples
#'
#' par(mfrow=c(2, 2))
#' obj1 <- ahead::fitforecast(AirPassengers)
#' obj2 <- ahead::fitforecast(AirPassengers, conformalize = TRUE)
#' plot(AirPassengers)
#' plot(obj1)
#' obj2$plot()
#' obj2$plot("simulations")
#'
fitforecast <- function(y,
                        h = NULL,
                        pct_train = 0.9,
                        pct_calibration = 0.5,
                        method = c("thetaf",
                                   "arima",
                                   "ets",
                                   "te",
                                   "tbats",
                                   "tslm",
                                   "dynrmf",
                                   "ridge2f",
                                   "naive",
                                   "snaive"),
                        level = 95,
                        B = 1000L,
                        seed = 17223L,
                        graph = TRUE,
                        conformalize = FALSE,
                        type_calibration = c("splitconformal",
                                             "cv1",
                                             "loocv"),
                        gap = 3L,
                        agg = c("mean", "median"),
                        vol = c("constant", "garch"),
                        type_sim = c("kde", "surrogate", "bootstrap"),
                        ...)
{
  method <- match.arg(method)
  agg <- match.arg(agg)
  vol <- match.arg(vol)
  type_sim <- match.arg(type_sim)
  type_calibration <- match.arg(type_calibration)

  if (is.null(h))
  {
    splitted_y_obj <- splitts(y, split_prob = pct_train)
    y_train <- splitted_y_obj$training
    y_test <- splitted_y_obj$testing
    h_test <- length(y_test)
  } else {
    y_train <- y
    h_test <- h
  }

  tspx <- tsp(y_train)
  start_preds <- tspx[2] + 1 / tspx[3]

  if (conformalize == TRUE)
  {
    splitted_y_train_obj <- splitts(y_train,
                                    split_prob = pct_calibration)
    y_train_calibration <- splitted_y_train_obj$training
    y_val_calibration <- splitted_y_train_obj$testing
    h_calibration <- length(y_val_calibration)
  }

  if (method %in% c("arima", "ets",
                    "tbats", "tslm"))
  {
    method <- match.arg(method)

    fcast_method <- switch(
      method,
      arima = forecast::auto.arima,
      ets = forecast::ets,
      tbats = forecast::tbats,
      tslm = function (y)
        tslm(y ~ trend + season)
    )

    fcast_func <- function(y)
    {
      return(do.call(what = fcast_method, args = list(y = y)))
    }

    if (conformalize == TRUE &&
        type_calibration == "splitconformal")
    {
      obj_calibration <-
        forecast::forecast(fcast_func(y = y_train_calibration),
                           h = h_calibration,
                           level = level)
      obj <- forecast::forecast(fcast_func(y = y_val_calibration),
                                h = h_test,
                                level = level)
    } else {
      obj <- forecast::forecast(fcast_func(y = y_train),
                                h = h_test,
                                level = level)
    }

  }

  if (method %in% c("dynrmf", "ridge2f",
                    "naive", "snaive",
                    "thetaf", "te"))
  {
    method <- match.arg(method)

    fcast_func <- switch(
      method,
      dynrmf = ahead::dynrmf,
      ridge2f = ahead::ridge2f,
      naive = forecast::naive,
      snaive = forecast::snaive,
      thetaf = forecast::thetaf,
      te = function (y, h, level, ...) {ahead::eatf(y, h, level, method = "EAT",
                                    weights = c(0.5, 0, 0.5), ...)}
    )

    if (conformalize == TRUE &&
        type_calibration == "splitconformal")
    {
      obj_calibration <- do.call(fcast_func,
                                 list(y = y_train_calibration,
                                      h = h_calibration,
                                      level = level))
      obj <- do.call(fcast_func,
                     list(
                       y = y_val_calibration,
                       h = h_test,
                       level = level,
                       ...
                     ))
    } else {
      obj <- do.call(fcast_func,
                     list(
                       y = y_train,
                       h = h_test,
                       level = level,
                       ...
                     ))
    }
  }

  out <- obj
  out$test <-
    stats::Box.test(out$residuals, lag = 1, type = "Ljung-Box") # examining the null hypothesis of independence in a given time series

  if (conformalize == TRUE)
  {
    stopifnot(!is.null(B) && B > 1)

    if (type_calibration == "splitconformal")
    {
      y_pred_calibration <- obj_calibration$mean
      calibrated_raw_residuals <-
        as.numeric(y_val_calibration) - as.numeric(y_pred_calibration)
    }

    if (type_calibration == "cv1")
    {
      n_train_calibration <- length(y_train_calibration)
      n_val_calibration <- length(y_val_calibration)
      y_pred_calibration <-
        calibrated_raw_residuals <- rep(0, n_val_calibration)
      obj_calibration <- do.call(fcast_func,
                                 list(y = y_train_calibration,
                                      h = 1L,
                                      level = level))
      pb <- utils::txtProgressBar(min = 0,
                                  max = n_val_calibration,
                                  style = 3)
      for (i in seq_along(n_val_calibration))
      {
        y_pred_calibration[i] <- as.numeric(obj_calibration$mean)
        calibrated_raw_residuals[i] <-
          y_val_calibration[i] - y_pred_calibration[i]
        obj_calibration <- do.call(fcast_func,
                                   list(
                                     y = c(y_train_calibration,
                                           y_val_calibration[seq_along(i)]),
                                     h = 1L,
                                     level = level
                                   ))
        utils::setTxtProgressBar(pb, i)
      }
      close(pb)
    }

    if (type_calibration == "loocv")
    {
      n_train_calibration <- length(y_train_calibration)
      n_val_calibration <- length(y_val_calibration)
      y_pred_calibration <-
        calibrated_raw_residuals <-
        matrix(0, nrow = n_val_calibration,
               ncol = n_val_calibration)
      pb <- utils::txtProgressBar(min = 0,
                                  max = n_val_calibration,
                                  style = 3)
      for (i in seq_along(n_val_calibration))
      {
        obj_calibration <- do.call(fcast_func,
                                   list(
                                     y = y_train_calibration[-(i:(i + gap))],
                                     h = h_calibration,
                                     level = level
                                   ))
        y_pred_calibration[i, ] <- as.numeric(obj_calibration$mean)
        calibrated_raw_residuals[i, ] <-
          y_val_calibration - y_pred_calibration[i, ]
        utils::setTxtProgressBar(pb, i)
      }
      close(pb)

    }

    if (type_calibration == "loocv")
    {
      y_pred_calibration <- switch(
        agg,
        median = apply(y_pred_calibration, 1, median),
        mean = apply(y_pred_calibration, 1, mean)
      )
      calibrated_raw_residuals <-
        rgaussiandens(as.vector(calibrated_raw_residuals),
                      n = h_calibration,
                      seed = seed)
    }

    matrix_preds_calibration <- replicate(B, y_pred_calibration)
    matrix_preds <- replicate(B, obj$mean)

    if (type_sim == "kde") {
      if (vol == "constant")
      {
        scaled_calibrated_residuals <- base::scale(calibrated_raw_residuals,
                                                   center = TRUE,
                                                   scale = TRUE)
        simulated_scaled_calibrated_residuals <-
          rgaussiandens(
            scaled_calibrated_residuals,
            n = h_test,
            p = B,
            seed = seed
          )
        sd_calibrated_residuals <- sd(calibrated_raw_residuals)
      } else {
        stopifnot(vol == "garch")
        fit_res <-
          suppressWarnings(fit_garch(calibrated_raw_residuals))
        omega <- fit_res@fit$coef["omega"]
        alpha <- fit_res@fit$coef["alpha1"]
        beta <- fit_res@fit$coef["beta1"]
        sigmat <- rep(0, length(calibrated_raw_residuals))
        sigmat[1] <- sqrt(mean(calibrated_raw_residuals ^ 2))
        for (i in 2:length(calibrated_raw_residuals))
          sigmat[i] <-
          sqrt(omega + alpha * sigmat[i - 1] ^ 2 + beta * calibrated_raw_residuals[i -
                                                                                     1] ^ 2)
        scaled_calibrated_residuals <-
          base::scale(calibrated_raw_residuals,
                      center = TRUE,
                      scale = FALSE) / sigmat
        simulated_scaled_calibrated_residuals <-
          rgaussiandens(
            scaled_calibrated_residuals,
            n = h_test,
            p = B,
            seed = seed
          )
        sd_calibrated_residuals <-
          matrix(0, nrow = h_test, ncol = B)
        sd_calibrated_residuals[1,] <-
          sqrt(omega + alpha * rev(sigmat)[1] ^ 2 + beta * rev(calibrated_raw_residuals)[1] ^
                 2)
        for (i in 2:h_test)
          sd_calibrated_residuals[i,] <-
          sqrt(
            omega + alpha * sd_calibrated_residuals[(i - 1),] ^ 2 + beta * simulated_scaled_calibrated_residuals[(i -
                                                                                                                    1),] ^ 2
          )
      }
    }

    if (type_sim == "surrogate") {
      if (vol == "constant")
      {
        scaled_calibrated_residuals <- base::scale(calibrated_raw_residuals,
                                                   center = TRUE,
                                                   scale = TRUE)
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
          tseries::surrogate(scaled_calibrated_residuals,
                                   ns =
                                     B)[seq_along(h_test), ]

        sd_calibrated_residuals <- sd(calibrated_raw_residuals)
      } else {
        stopifnot(vol == "garch")
        fit_res <-
          suppressWarnings(fit_garch(calibrated_raw_residuals))
        omega <- fit_res@fit$coef["omega"]
        alpha <- fit_res@fit$coef["alpha1"]
        beta <- fit_res@fit$coef["beta1"]
        sigmat <- rep(0, length(calibrated_raw_residuals))
        sigmat[1] <- sqrt(mean(calibrated_raw_residuals ^ 2))
        for (i in 2:length(calibrated_raw_residuals))
          sigmat[i] <-
          sqrt(omega + alpha * sigmat[i - 1] ^ 2 + beta * calibrated_raw_residuals[i -
                                                                                     1] ^ 2)
        scaled_calibrated_residuals <-
          base::scale(calibrated_raw_residuals,
                      center = TRUE,
                      scale = FALSE) / sigmat
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
          tseries::surrogate(scaled_calibrated_residuals,
                                   ns =
                                     B)[seq_along(h_test), ]
        sd_calibrated_residuals <-
          matrix(0, nrow = h_test, ncol = B)
        sd_calibrated_residuals[1,] <-
          sqrt(omega + alpha * rev(sigmat)[1] ^ 2 + beta * rev(calibrated_raw_residuals)[1] ^
                 2)
        for (i in 2:h_test)
          sd_calibrated_residuals[i,] <-
          sqrt(
            omega + alpha * sd_calibrated_residuals[(i - 1),] ^ 2 + beta * simulated_scaled_calibrated_residuals[(i -
                                                                                                                    1),] ^ 2
          )
      }
    }

    if (type_sim == "bootstrap") {
      freq_calibrated_raw_residuals <- frequency(calibrated_raw_residuals)
      if (length(calibrated_raw_residuals) <= 2 * freq_calibrated_raw_residuals)
        freq_calibrated_raw_residuals <- 1L
      block_size <-
        ifelse(
          freq_calibrated_raw_residuals > 1,
          2 * freq_calibrated_raw_residuals,
          min(8, floor(
            length(calibrated_raw_residuals) / 2
          ))
        )
      block_size <-
        floor(min(
          max(3L, block_size),
          length(calibrated_raw_residuals) - 1L
        ))

      if (vol == "constant")
      {
        scaled_calibrated_residuals <- base::scale(calibrated_raw_residuals,
                                                   center = TRUE,
                                                   scale = TRUE)
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
          tseries::tsbootstrap(
            scaled_calibrated_residuals,
            nb =
              B,
            b = floor(block_size),
            type =
              "block"
          )[seq_along(h_test), ]
        sd_calibrated_residuals <- sd(calibrated_raw_residuals)
      } else {
        stopifnot(vol == "garch")
        fit_res <-
          suppressWarnings(fit_garch(calibrated_raw_residuals))
        omega <- fit_res@fit$coef["omega"]
        alpha <- fit_res@fit$coef["alpha1"]
        beta <- fit_res@fit$coef["beta1"]
        sigmat <- rep(0, length(calibrated_raw_residuals))
        sigmat[1] <- sqrt(mean(calibrated_raw_residuals ^ 2))
        for (i in 2:length(calibrated_raw_residuals))
          sigmat[i] <-
          sqrt(omega + alpha * sigmat[i - 1] ^ 2 + beta * calibrated_raw_residuals[i -
                                                                                     1] ^ 2)
        scaled_calibrated_residuals <-
          base::scale(calibrated_raw_residuals,
                      center = TRUE,
                      scale = FALSE) / sigmat
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
          tseries::tsbootstrap(
            scaled_calibrated_residuals,
            nb =
              B,
            b = block_size,
            type =
              "block"
          )[seq_along(h_test), ]
        sd_calibrated_residuals <-
          matrix(0, nrow = h_test, ncol = B)
        sd_calibrated_residuals[1,] <-
          sqrt(omega + alpha * rev(sigmat)[1] ^ 2 + beta * rev(calibrated_raw_residuals)[1] ^
                 2)
        for (i in 2:h_test)
          sd_calibrated_residuals[i,] <-
          sqrt(
            omega + alpha * sd_calibrated_residuals[(i - 1),] ^ 2 + beta * simulated_scaled_calibrated_residuals[(i -
                                                                                                                    1),] ^ 2
          )
      }
    }

    sims <-
      matrix_preds + sd_calibrated_residuals * simulated_scaled_calibrated_residuals

    sims <- ts(sims,
               start = start_preds,
               frequency = frequency(y_test))
    preds_lower <-
      apply(sims, 1, function(x)
        quantile(x, probs = (1 - level / 100) / 2))
    preds_upper <-
      apply(sims, 1, function(x)
        quantile(x, probs = 1 - (1 - level / 100) / 2))
    out$model <- list(obj = obj,
                      obj_calibration = obj_calibration)
    out$method <- paste0("\nconformalized ", out$method)
    out$mean <- ts(switch(
      agg,
      median = apply(sims, 1, median),
      mean = apply(sims, 1, mean)
    ),
    start = start(y_test),
    frequency = frequency(y_test))
    out$lower <- ts(preds_lower,
                    start = start(y_test),
                    frequency = frequency(y_test))
    out$upper <- ts(preds_upper,
                    start = start(y_test),
                    frequency = frequency(y_test))
    out$x <- y_train
    out$residuals <- ts(
      calibrated_raw_residuals,
      start = start(y_val_calibration),
      frequency = frequency(y_val_calibration)
    ) # /!\ not the same residuals, beware
    out$fitted <- ts(
      y_pred_calibration,
      start = start(y_val_calibration),
      frequency = frequency(y_val_calibration)
    ) # /!\ not the same fitted, beware
    out$sims <- sims
    out$test <-
      stats::Box.test(out$residuals, type = "Ljung-Box") # examining the null hypothesis of independence in a given time series
    out <- structure(out, class = "forecast")
  }

  out$coverage <- 100 * mean((y_test >= out$lower) * (out$upper >= y_test))
  out$interval_length <- mean(out$upper - out$lower)
  out$winkler_score <- winkler_score2(out, y_test, level = level)
  accuracy_results <- forecast::accuracy(out, y_test)
  out$ME <- accuracy_results[2, "ME"]
  out$RMSE <- accuracy_results[2, "RMSE"]
  out$MAE <- accuracy_results[2, "MAE"]
  out$MPE <- accuracy_results[2, "MPE"]
  out$MAPE <- accuracy_results[2, "MAPE"]
  out$MASE <- accuracy_results[2, "MASE"]

  plots <-
    function(out,
             type = c("forecast", "simulations", "residuals"))
    {
      if (missing(type))
      {
        type <- "forecast"
      } else {
        type <- match.arg(type)
      }

      if (type == "forecast")
      {
        plot(out, lwd = 2)
        lines(y_test,
              col = "red",
              lwd = 2,
              lty = 2)
      }

      if (type == "simulations")
      {
        matplot(
          x = as.numeric(time(y)),
          y = y,
          type = "l",
          lwd = 2,
          main = paste0(B, " predictive simulations for ", out$method, "\n"),
          xlab = "",
          ylab = "",
          ylim = range(c(y, out$sims))
        )
        matlines(x = as.numeric(time(y_test)),
                 y = out$sims,
                 lwd = 2)
      }

      if (type == "residuals")
      {
        forecast::checkresiduals(out)
      }
    }

  if (graph == TRUE)
  {
    out$plot <- function (type)
      plots(out, type)
  }

  return(out)
}
