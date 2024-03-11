# 0 - interface function for dyrm -----------------------------------------

#' Dynamic regression model
#'
#' Adapted from forecast::nnetar, with alternative fitting functions (see examples)
#'
#' @param y A numeric vector or time series of class \code{ts}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param fit_func Fitting function (Statistical/ML model). Default is Ridge regression.
#' @param predict_func Prediction function (Statistical/ML model)
#' @param fit_params a list of additional parameters for the fitting function \code{fit_func} (see examples)
#' @param type_pi Type of prediction interval (currently "gaussian", ETS: "E", Arima: "A" or Theta: "T")
#' @param xreg_fit Optionally, a vector or matrix of external regressors, which
#' must have the same number of rows as y. Must be numeric.
#' @param xreg_predict Future values of external regressor variables.
#' @param ... additional parameters
#'
#' @return a list; an object of class \code{forecast}.
#'
#' The function \code{summary} is used to obtain and print a summary of the
#' results.
#'
#' The generic accessor functions \code{fitted.values} and \code{residuals}
#' extract useful features.
#'
#' @author T. Moudiki
#'
#' @export
#'
#' @references
#'
#' Hyndman, R. J., & Athanasopoulos, G. (2018). Forecasting:
#' principles and practice. OTexts.
#'
#' Hyndman R, Athanasopoulos G, Bergmeir C, Caceres G, Chhay L,
#' O'Hara-Wild M, Petropoulos F, Razbash S, Wang E, Yasmeen F (2021).
#' forecast: Forecasting functions for time series and linear models. R
#' package version 8.14, <URL: https://pkg.robjhyndman.com/forecast/>. \cr
#'
#'
#' @examples
#'
#' # Example 0: with Ridge regression
#'
#' par(mfrow=c(3, 2))
#' plot(dynrmf(USAccDeaths, h=20, level=95))
#' plot(dynrmf(AirPassengers, h=20, level=95))
#' plot(dynrmf(lynx, h=20, level=95))
#' plot(dynrmf(WWWusage, h=20, level=95))
#' plot(dynrmf(Nile, h=20, level=95))
#' plot(dynrmf(fdeaths, h=20, level=95))
#'
#'
#' # Example 1: with Random Forest
#'
#' \dontrun{
#'
#' require(randomForest)
#'
#' par(mfrow=c(3, 2))
#' plot(dynrmf(USAccDeaths, h=20, level=95, fit_func = randomForest::randomForest,
#'      fit_params = list(ntree = 50), predict_func = predict))
#' plot(dynrmf(AirPassengers, h=20, level=95, fit_func = randomForest::randomForest,
#'      fit_params = list(ntree = 50), predict_func = predict))
#' plot(dynrmf(lynx, h=20, level=95, fit_func = randomForest::randomForest,
#'      fit_params = list(ntree = 50), predict_func = predict))
#' plot(dynrmf(WWWusage, h=20, level=95, fit_func = randomForest::randomForest,
#'      fit_params = list(ntree = 50), predict_func = predict))
#' plot(dynrmf(Nile, h=20, level=95, fit_func = randomForest::randomForest,
#'      fit_params = list(ntree = 50), predict_func = predict))
#' plot(dynrmf(fdeaths, h=20, level=95, fit_func = randomForest::randomForest,
#'      fit_params = list(ntree = 50), predict_func = predict))
#'}
#'
#' # Example 2: with SVM
#'
#'\dontrun{
#'
#' require(e1071)
#'
#' par(mfrow=c(2, 2))
#' plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
#' fit_params = list(kernel = "linear"), predict_func = predict))
#' plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
#' fit_params = list(kernel = "polynomial"), predict_func = predict))
#' plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
#' fit_params = list(kernel = "radial"), predict_func = predict))
#' plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
#' fit_params = list(kernel = "sigmoid"), predict_func = predict))
#'
#'}
#'
#'
dynrmf <- function(y, h = 5,
                   level = 95,
                   fit_func = ahead::ridge,
                   predict_func = predict,
                   fit_params = NULL,
                   type_pi = c("gaussian", "E", "A", "T"),
                   xreg_fit = NULL,
                   xreg_predict = NULL,
                   ...)
{
  if(is_package_available("forecast") == FALSE)
  {
    try_forecast <- try(utils::install.packages("forecast",
                                                repos = c("https://cloud.r-project.org",
                                                          "https://cran.rstudio.com")),
    silent = TRUE)
    if(inherits(try_forecast, "try-error"))
    {
      try(utils::install.packages("remotes",
                                  repos = c("https://cloud.r-project.org",
                                            "https://cran.rstudio.com")))
      try(remotes::install_github('robjhyndman/forecast',
                                  dependencies=TRUE),
          silent=TRUE)
    }
  }

  stopifnot(length(level) == 1)

  type_pi <- match.arg(type_pi)

  obj <- dynrm_fit(y,
                   fit_func = fit_func,
                   predict_func = predict_func,
                   fit_params = fit_params,
                   xreg = xreg_fit,
                   ...)

  # in predict: return(structure(out, class = "forecast"))
  return(dynrm_predict(obj, h=h, level=level, type_pi=type_pi,
                       xreg = xreg_predict,
                       ...))
}
dynrmf <- compiler::cmpfun(dynrmf)



# 1 - dynrm fitting ---------------------------------------------------------------

# Fit a dynrm (adapted from forecast::nnetar)
dynrm_fit <- function(y,
                      p,
                      P = 1,
                      xreg = NULL,
                      fit_func = ahead::ridge,
                      predict_func = predict,
                      fit_params = NULL,
                      lambda = NULL, alpha = NULL,
                      scale_inputs = TRUE,
                      x = y,
                      centers = NULL,
                      seed = 123,
                      ...) {

  yname <- deparse(substitute(y))

  # transform data
  if (!is.null(lambda)) {
    xx <- forecast::BoxCox(x, lambda)
    lambda <- attr(xx, "lambda")
  } else {
    xx <- x
  }

  # scale series
  scalex <- NULL
  if (scale_inputs) {
    tmpx <- scale_ahead(xx, center = TRUE, scale = TRUE)
    scalex <- list(
      center = attr(tmpx, "scaled:center"),
      scale = attr(tmpx, "scaled:scale")
    )
    xx <- scale_ahead(xx, center = scalex$center,
                      scale = scalex$scale)
    xx <- xx[, 1]
  }


  # check xreg class & dim
  xxreg <- NULL
  scalexreg <- NULL
  if (!is.null(xreg)) {
    xxreg <- xreg <- as.matrix(xreg)

    if (length(x) != NROW(xreg)) {
      stop("Number of rows in xreg does not match series length")
    }

    # scale xreg
    if (scale_inputs) {
      tmpx <- base::scale(xxreg, center = TRUE, scale = TRUE)
      scalexreg <- list(
        center = attr(tmpx, "scaled:center"),
        scale = attr(tmpx, "scaled:scale")
      )

      xxreg <- scale(xxreg,
                     center = scalexreg$center,
                     scale = scalexreg$scale)
    }

  }


  # Set up lagged matrix
  n <- length(xx)
  xx <- as.ts(xx)
  m <- max(round(frequency(xx)), 1L)

  if (m == 1) {

    if (missing(p)) {
      p <- max(length(stats::ar(xx)$ar), 1)
    }

    if (p >= n) {
      warning("Reducing number of lagged inputs due to short series")
      p <- n - 1
    }

    lags <- 1:p

    if (P > 1) {
      warning("Non-seasonal data, ignoring seasonal lags")
    }

    P <- 0

  } else { # if m != 1

    if (missing(p)) {
      if (n > 2 * m) {
        x.sa <- forecast::seasadj(forecast::mstl(forecast::na.interp(xx)))
      } else {
        x.sa <- forecast::na.interp(xx)
      }

      p <- max(length(stats::ar(x.sa)$ar), 1)
    }

    if (p >= n) {
      warning("Reducing number of lagged inputs due to short series")
      p <- n - 1
    }

    if (P > 0 && n >= m * P + 2) {
      lags <- sort(unique(c(1:p, m * (1:P))))
    } else {
      lags <- 1:p
      if (P > 0) {
        warning("Series too short for seasonal lags")
        P <- 0
      }
    }

  }

  maxlag <- max(lags)
  nlag <- length(lags)
  y <- xx[-(1:maxlag)]
  lags.X <- matrix(NA_real_, ncol = nlag, nrow = n - maxlag)

  for (i in 1:nlag)
    lags.X[, i] <- xx[(maxlag - lags[i] + 1):(n - lags[i])]

  if (is.null(alpha) == FALSE)
  {
    stop("not used yet")
    # n_alphas <- nrow(lags.X)
    # alphas <- alpha*((1-alpha)^((n_alphas + p + P - 2):0))
    # # Add xreg into lagged matrix
    # lags.X <- cbind(lags.X*exp(embedc(alphas, ncol(lags.X)-1)),
    #                 xxreg[-(1:maxlag), ])
  } else {
    # Add xreg into lagged matrix
    lags.X <- cbind(lags.X, xxreg[-(1:maxlag), ])
  }

  # cat("lags.X: ", "\n")
  # print(head(lags.X))
  # print(tail(lags.X))
  # cat("\n")

  # Remove missing values if present
  j <- complete.cases(lags.X, y)
  X_j <- as.matrix(lags.X[j, , drop = FALSE])

  # if (!is.null(centers))
  # {
  #   # cat("X_j", "\n")
  #   # print(head(X_j))
  #   # cat("\n")
  #   if(ncol(X_j) > 1)
  #   {
  #     # print(dim(X_j))
  #      print("here")
  #     fit_cclust <- cclust::cclust(X_j, centers = centers,
  #                                  iter.max = 20)
  #     X_j <- cbind(X_j, fit_cclust$cluster)
  #   } else {
  #      print("here2")
  #     # cat("X_j", "\n")
  #     # print(head(X_j))
  #     # cat("\n")
  #     (fit_cclust <- cclust::cclust(cbind(rep(1L, nrow(X_j)), X_j),
  #                                   centers = centers,
  #                                   iter.max = 20))
  #     X_j <- cbind(X_j, fit_cclust$cluster)
  #     # cat("X_j", "\n")
  #     # print(head(X_j))
  #     # cat("\n")
  #   }
  #
  # }

  ## Stop if there's no data to fit (e.g. due to NAs or NaNs)
  if (NROW(X_j) == 0) {
    stop("No data to fit (possibly due to NA or NaN)")
  }

  # cat("y: ", "\n")
  # print(head(y))
  # print(tail(y))
  # cat("\n")

  # fit
  set.seed(seed) # in case the algo in fit_func is randomized

  # cat("X_j", "\n")
  # print(X_j)
  # cat("\n")

  fit <- do.call(what = fit_func,
                 args = c(list(x = X_j,
                               y = y[j]), fit_params))

   # cat("fit", "\n")
   # print(fit)
   # cat("\n")

  # Return results
  out <- list()
  out$x <- as.ts(x)
  out$m <- m
  out$p <- p
  out$P <- P
  out$scalex <- scalex
  out$scalexreg <- scalexreg
  out$xreg <- xreg
  out$lambda <- lambda
  out$model <- fit
  out$dynrmgs <- list(...)


   # cat("fit", "\n")
   # print(fit)
   # cat("\n")


    #cat("head(X_j)", "\n")
    #print(head(X_j))
    #cat("\n")

  # Fitted values
  fits <- try(drop(predict_func(fit, newdata = X_j)),
              silent = TRUE)
  if (inherits(fits, "try-error") || is.null(fits))
  {
     fits <- try(drop(predict_func(fit,
                              newx = X_j)),
                 silent = TRUE)
     if (inherits(fits, "try-error") || is.null(fits))
     {
       fits <- try(drop(predict_func(fit, X_j)),
                   silent = TRUE)
       }
   }

   # cat("fits", "\n")
   # print(fits)
   # cat("\n")

  if (scale_inputs) {
      # cat("scalex$scale", "\n")
      # print(scalex$scale)
      # cat("\n")
      # cat("scalex$center", "\n")
      # print(scalex$center)
      # cat("\n")
    fits <- fits * scalex$scale + scalex$center
  }

  fits <- ts(fits, end = end(out$x),
             frequency = frequency(out$x))
  if (!is.null(lambda)) {
    fits <- forecast::InvBoxCox(fits, lambda)
  }

  out$fitted <- fits
  out$residuals <- out$x - out$fitted
  out$lags <- lags
  out$predict_func <- predict_func
  out$series <- yname
  out$method <- paste("DynRM ", p, sep = "")
  out$sigma <- sd(out$residuals)
  # if (!is.null(centers))
  # {
  #   out$centers <- centers
  #   out$fit_cclust <- fit_cclust
  # }

  if (P > 0) {
    out$method <- paste(out$method, ",", P, sep = "")
  }

  if (P > 0) {
    out$method <- paste(out$method, "[", m, "]", sep = "")
  }

  out$call <- match.call()

  # return
  invisible(structure(out, class = c("dynrm")))
}


# 2 - dynrm prediction ----------------------------------------------------

# Predict from dynrm
dynrm_predict <- function(out,
                          h = ifelse(out$m > 1, 2 * out$m, 10),
                          level = 95,
                          fan = FALSE,
                          xreg = NULL,
                          lambda = out$lambda,
                          type_pi = c("gaussian", "E", "A", "T"),
                          ...) {

  # cat("\n")
  # print("in dynrm_predict")

  # print("out$x")
  # print(out$x)
  # print("\n")

  tspx <- tsp(out$x)

  type_pi <- match.arg(type_pi)

  # print("tspx")
  # print(tspx)
  # print("\n")

  if (fan) {
    level <- seq(51, 99, by = 3)
  } else {
    if (min(level) > 0 && max(level) < 1) {
      level <- 100 * level
    } else if (min(level) < 0 || max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }


  # Check if xreg was used in fitted model
  if (is.null(out$xreg)) {
    if (!is.null(xreg)) {
      warning("External regressors were not used in fitted model, xreg will be ignored")
    }
    xreg <- NULL
  }
  else {
    if (is.null(xreg)) {
      stop("No external regressors provided")
    }
    xreg <- as.matrix(xreg)
    if (NCOL(xreg) != NCOL(out$xreg)) {
      stop("Number of external regressors does not match fitted model")
    }
    if (!identical(colnames(xreg), colnames(out$xreg))) {
      warning(
        "xreg contains different column names from the xreg used in training. Please check that the regressors are in the same order."
      )
    }
    h <- NROW(xreg)
  }


  fcast <- numeric(h)
  xx <- out$x
  xxreg <- xreg
  if (!is.null(lambda)) {
    xx <- forecast::BoxCox(xx, lambda)
    lambda <- attr(xx, "lambda")
  }


  # Check and apply scaling of fitted model
  if (!is.null(out$scalex)) {
    # scale_ahead is in utils.R
    xx <- scale_ahead(xx,
                      center = out$scalex$center,
                      scale = out$scalex$scale)
    if (!is.null(xreg)) {
      xxreg <- base::scale(xreg,
                           center = out$scalexreg$center,
                           scale = out$scalexreg$scale)
    }
  }


  # Get lags used in fitted model
  lags <- out$lags
  maxlag <- max(lags)
  flag <- rev(tail(xx, n = maxlag))


  if (is.null(out$centers))
  {
    # Iterative 1-step forecast
    for (i in 1:h)
    {
      newdata <- c(flag[lags], xxreg[i, ])


      newdata_ <- as.matrix(rbind(newdata,
                                  rnorm(length(newdata)))) # do not change (ever)

      # cat("newdata", "\n")
      # print(newdata_)
      # cat("\n")

      preds <- try(as.numeric(out$predict_func(out$model,
                                               newdata = newdata_)[1]), # do not change (ever)
                   silent = TRUE)
      if (inherits(preds, "try-error"))
      {
        preds <- try(as.numeric(out$predict_func(out$model,
                                                 newx = newdata_)[1]), # do not change (ever)
                     silent = TRUE)
        if (inherits(preds, "try-error"))
        {
          preds <- NA
        }
      }

      # cat("preds", "\n")
      # print(preds)
      # cat("\n")

      fcast[i] <- preds

      flag <- c(fcast[i], flag[-maxlag])
    }

  } else { # centers not NULL # TODO

    stop("clustering not implemented yet")
  }

  # Re-scale point forecasts
  if (!is.null(out$scalex)) {
    #fcast <- fcast * out$scalex$scale + out$scalex$center
    fcast <- fcast * out$scalex$scale + out$scalex$center
  }


  # Add ts properties
  fcast <- ts(fcast, start = tspx[2] + 1 / tspx[3],
              frequency = tspx[3])


  # Back-transform point forecasts
  if (!is.null(lambda)) {
    fcast <- forecast::InvBoxCox(fcast, lambda)
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
    rep_qt_sd <- ts(rep(qt_sd, h), start = tspx[2] + 1 / tspx[3],
                    frequency = tspx[3])

    resid_fcast <- list()
    resid_fcast$mean <- 0
    resid_fcast$lower <- -rep_qt_sd
    resid_fcast$upper <- rep_qt_sd
  }

  out$mean <-  drop(fcast  + resid_fcast$mean)
  out$lower <- drop(fcast + resid_fcast$lower)
  out$upper <- drop(fcast + resid_fcast$upper)
  out$level <- level

  return(structure(out, class = "forecast"))
}
