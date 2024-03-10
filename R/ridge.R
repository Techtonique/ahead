#' #' Forecasts from ridge regression models
#' #'
#' #' @param y A multivariate (for now) time series of class \code{ts} (preferred) or a \code{matrix}
#' #' @param h Forecasting horizon
#' #' @param level Confidence level for prediction intervals
#' #' @param lags Number of lags to include in the model
#' #' @param nb_hidden Number of hidden units in the neural network
#' #' @param fit_method Method to fit the model. Either \code{ridge} for elastic net or \code{mgaussian} for multivariate Gaussian elastic net
#' #' @param nodes_sim Method to generate the nodes for the simulation of the multivariate . Either \code{sobol}, \code{halton} or \code{unif}
#' #' @param activ Activation function for the hidden layer. Either \code{relu}, \code{sigmoid}, \code{tanh}, \code{leakyrelu}, \code{elu} or \code{linear}
#' #' @param hidden_layer_bias Whether to include a bias term in the hidden layer
#' #' @param a Hyperparameter for activation function "leakyrelu", "elu"
#' #' @param lambda Hyperparameter for elastic net (regularization parameter)
#' #' @param alpha Hyperparameter for elastic net (compromise between ridge and lasso)
#' #' @param seed Seed for reproducibility
#' #' @param type_ci Type of confidence interval. Either \code{none}, \code{parametric} or \code{nonparametric}
#' #' @param type_forecast Type of forecast. Either \code{recursive} or \code{direct}
#' #'
#' #' @return an object of class \code{mforecast}; a list containing the following elements:
#' #' \item{mean}{Point forecasts for the time series} \item{resid}{Residuals from the fitted model}
#' #' \item{model}{A list containing information about the fitted model}
#' #' \item{method}{The name of the forecasting method as a character string}
#' #' \item{x}{The original time series}
#' #'
#' #' @export
#' #'
#' #' @examples
#' #'
#' #'
#' ridgef <- function(y, h = 5,
#'                    level = 95,
#'                    lags = 1,
#'                    nb_hidden = 5,
#'                    fit_method = c("ridge", "mgaussian"),
#'                    nodes_sim = c("sobol", "halton", "unif"),
#'                    activ = c("relu", "sigmoid", "tanh",
#'                              "leakyrelu", "elu", "linear"),
#'                    hidden_layer_bias = FALSE,
#'                    a = 0.01,
#'                    lambda = 0.1,
#'                    alpha = 0.5,
#'                    seed = 1,
#'                    type_pi = "none",
#'                    type_forecast = c("recursive", "direct"))
#' {
#'   if (!is.ts(y))
#'   {
#'     x <- ts(y)
#'   } else {
#'     x <- y
#'   }
#'
#'   freq_x <- frequency(x)
#'   start_fits <- start(x)
#'   start_preds <- tsp(x)[2] + 1/freq_x
#'   n_series <- ncol(x)
#'
#'   fit_method <- match.arg(fit_method)
#'   nodes_sim <- match.arg(nodes_sim)
#'   activ <- match.arg(activ)
#'   type_forecast <- match.arg(type_forecast)
#'
#'   # Fitting a regularized regression model to multiple time series
#'   fit_obj <- fit_ridge_mts(
#'     x,
#'     lags = lags,
#'     nb_hidden = nb_hidden,
#'     fit_method = fit_method,
#'     nodes_sim = nodes_sim,
#'     activ = activ,
#'     hidden_layer_bias = hidden_layer_bias,
#'     a = a,
#'     lambda = lambda,
#'     alpha = alpha,
#'     seed = seed
#'   )
#'
#'   # Forecast from fit_obj
#'
#'   preds <- ts(data = fcast_obj_mts(fit_obj,
#'     h = h,
#'     type_pi = type_pi,
#'     type_forecast = type_forecast,
#'     level = level),
#'     start = start_preds, frequency = freq_x)
#'
#'   resids <- ts(data = fit_obj$resid,
#'                start = start_preds, frequency = freq_x)
#'
#'   # Forecast from fit_obj
#'   forecasts <- lapply(1:n_series,
#'                       FUN = function (i) {structure(ts(preds[,i],
#'                       start = start_preds,
#'                       frequency = freq_x), class = "forecast")})
#'
#'   names(forecasts) <- colnames(x)
#'
#'   out <- list(mean = preds,
#'               residuals = resids,
#'               method = fit_method,
#'               model = fit_obj,
#'               x = x,
#'               forecast = forecasts)
#'
#'   return(structure(out, class = "mforecast"))
#' }
#' ridgef <- compiler::cmpfun(ridgef)
