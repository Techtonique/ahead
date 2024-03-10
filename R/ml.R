#'
#' # y <- AirPassengers;
#' # obj <- model_fit(y, lags=2, obj = ahead::ridge);
#' # (res <- model_predict(obj, h=20L));
#' # plot(res, type = 'l')
#'
#' # y <- Nile;
#' # obj <- model_fit(y, lags=2, obj = ahead::ridge);
#' # (res <- model_predict(obj, h=20L));
#' # plot(res, type = 'l')
#' #
#' # y <- USAccDeaths;
#' # obj <- model_fit(y, lags=2, obj = ahead::ridge);
#' # (res <- model_predict(obj, h=20L));
#' # plot(res, type = 'l')
#' #
#' # y <- mdeaths;
#' # obj <- model_fit(y, lags=2, obj = ahead::ridge);
#' # (res <- model_predict(obj, h=20L));
#' # plot(res, type = 'l')
#' #
#' # y <- austres
#' # obj <- model_fit(y, lags=2, obj = ahead::ridge);
#' # (res <- model_predict(obj, h=20L));
#' # plot(res, type = 'l')
#'
#'
#' #' Forecasting using Machine Leaning models
#' #'
#' #' @param y A numeric vector or time series of class \code{ts}
#' #' @param h Forecasting horizon
#' #' @param level Confidence level for prediction intervals
#' #' @param lags Number of lags of the input time series considered in the regression
#' #' @param fit_func Fitting function (Statistical/ML model). Default is Ridge regression.
#' #' @param predict_func Prediction function (Statistical/ML model)
#' #' @param type_pi Type of prediction interval
#' #' @param ... additional parameters passed to the fitting function \code{fit_func}
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' mlf <- function(y, h = 5,
#'                 level = 95,
#'                 lags = 1L,
#'                 fit_func = ahead::ridge,
#'                 predict_func = predict,
#'                 fit_params = NULL,
#'                 type_pi = c("gaussian", "E", "A", "T",
#'                             "kde", "surogate",
#'                             "bootstrap",
#'                             "blockbootsrap",
#'                             "movingblockbootsrap"),
#'                 ...)
#' {
#'
#'   type_pi <- match.arg(type_pi)
#'
#'   fit_obj <- model_fit(y, lags = lags,
#'                        obj = fit_func,
#'                        predict_func = predict_func,
#'                         ...)
#'   out <- model_predict(fit_obj, h = h,
#'                        predict_func = predict_func,
#'                        type_pi = type_pi,
#'                        level = level,
#'                        ...)
#'   out$x <- y
#'   out$method <- paste0(deparse(substitute(fit_func)),
#'                        paste0("(", out$lags, ")"))
#'   return(out)
#' }
#'
#' # 0 - inputs and parameters -----
#' # For uncertainty:
#' # Gaussian
#' # Independant bootstrap
#' # MB bootstrap
#' # circular block bootstrap
#' # rgaussiandens
#' # rsurr
#' # rboot
#' model_fit <- function(y,
#'                       lags = 1L,
#'                       obj = randomForest::randomForest,
#'                       predict_func = predict,
#'                       ...) {
#'   emb_y <- ahead::embedc(y, lags + 1L)
#'   p <- ncol(emb_y)
#'   if (lags == 1L) {
#'     x <- matrix(emb_y[,-p], ncol = 1)
#'   } else {
#'     x <- emb_y[,-p]
#'   }
#'   target <- ts(emb_y[, p], start=start(y),
#'                frequency=frequency(y))
#'   out <- list()
#'   out$lags <- lags
#'   out$x <- y
#'   scales <- scale(x)
#'   out$xm <- attributes(scales)$`scaled:center`
#'   out$xsd <- attributes(scales)$`scaled:scale`
#'   out$reg <- scales[,]
#'   out$ym <- mean(y)
#'   centered_y <- target - out$ym
#'   out$y <- centered_y
#'   debug_print(centered_y)
#'   debug_print(out$reg)
#'   out$obj <- obj(out$reg, centered_y, ...)
#'   out$fitted <- ts(predict_func(out$obj,
#'                                 out$reg),
#'                    start = start(y),
#'                    frequency = frequency(y)) + out$ym
#'   out$residuals <- out$y - out$fitted
#'   return(out)
#' }
#'
#' model_predict <- function(obj, h = 5,
#'                           level = 95,
#'                           predict_func = predict,
#'                           type_pi = c("gaussian", "E", "A", "T",
#'                                       "kde", "surogate",
#'                                       "bootstrap",
#'                                       "blockbootsrap",
#'                                       "movingblockbootsrap"),
#'                           ...) {
#'   type_pi <- match.arg(type_pi)
#'   lags <- obj$lags
#'   p <- ncol(obj$reg)
#'   n_y <- length(obj$y)
#'   target <- obj$y
#'   debug_print(target)
#'   tspy <- tsp(obj$x)
#'   start_preds <- tspy[2] + 1 / tspy[3]
#'   freq_preds <- frequency(obj$x)
#'   if (lags <= 1L) {
#'     newx <- matrix(target, ncol = 1)
#'   } else {
#'     newx <- cbind(obj$reg[, (p - lags + 2):p], target)
#'     colnames(newx) <- NULL
#'     newx <- sweep(newx, 2, obj$xm, "-")
#'     newx <- sweep(newx, 2, obj$xsd, "/")
#'   }
#'   debug_print(tail(newx))
#'
#'   target <- c(target, predict_func(obj$obj, newx, ...)[1] + obj$ym)
#'   debug_print(target)
#'   for (i in 2:h) {
#'     if (lags <= 1L) {
#'       newx <- matrix(target, ncol = 1)
#'     } else {
#'       newx <- cbind(newx[, (p - lags + 2):p], target)
#'       colnames(newx) <- NULL
#'     }
#'     newx <- sweep(newx, 2, obj$xm, "-")
#'     newx <- sweep(newx, 2, obj$xsd, "/")
#'     debug_print(tail(newx))
#'     target <- c(target, predict_func(obj$obj, newx, ...)[1] + obj$ym)
#'   }
#'   debug_print(target)
#'   out <- list()
#'   out$model <- obj
#'   out$x <- obj$x
#'   out$mean <- ts(target[(n_y + 1):(n_y + h)],
#'                  start = start_preds,
#'                  frequency = freq_preds)
#'   out$level <- level
#'   return(compute_uni_pi(out, h=h, type_pi=type_pi))
#' }
