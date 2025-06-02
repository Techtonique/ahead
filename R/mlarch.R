#' Conformalized Forecasting using Machine Learning models with ARCH effects
#'
#' @param y A numeric vector or time series of class \code{ts}
#' @param h Forecasting horizon
#' @param mean_model Function to fit the mean model (default: \code{forecast::auto.arima})
#' @param model_residuals Function to model the residuals (default: \code{forecast::thetaf})
#' @param fit_func Fitting function for the variance model (default: \code{ahead::ridge})
#' @param predict_func Prediction function for the variance model (default: \code{predict})
#' @param type_pi Type of prediction interval ("kde", "surrogate", or "bootstrap") for volatility modeling
#' @param type_sim_conformalize Type of simulation for conformalization of standardized residuals ("block-bootstrap", "surrogate", "kde", "bootstrap", or "fitdistr")
#' @param ml_method Machine learning method to use (if NULL, uses default methods)
#' @param level Confidence level for prediction intervals
#' @param B Number of bootstrap replications or simulations
#'
#' @return A forecast object containing predictions and prediction intervals
#' @export
#'
#' @examples
#' # Example usage will be added
mlarchf <- function(y, h=10L,
                    mean_model=forecast::auto.arima,
                    model_residuals=forecast::thetaf,
                    fit_func=ahead::ridge,
                    predict_func=predict,
                    type_pi = c("kde", "surrogate", "bootstrap"),
                    type_sim_conformalize = c("block-bootstrap", "surrogate", "kde", "bootstrap", "fitdistr"),
                    ml_method=NULL,
                    level=95,
                    B=250L)
{
  if (!is.ts(y))
  {
    y <- ts(y)
  }
  x <- y
  last_price <- as.numeric(tail(y, 1))
  freq_x <- frequency(y)
  start_preds <- tsp(y)[2] + 1/freq_x
  type_pi <- match.arg(type_pi)
  type_sim_conformalize <- match.arg(type_sim_conformalize)
  fit_mean <- ahead::genericforecast(mean_model, y=y, h=h)
  
  if (!is.null(ml_method))
  {
      fit_func <- function(x, y, method = ml_method, ...)
    {
      df <- data.frame(y=y, as.matrix(x)) 
      colnames(df) <- c("y", paste0("X", 1:ncol(x)))
      caret::train(y ~ ., data=df,
                  method = ml_method,
                  trControl=caret::trainControl(method = "none"))
    }

    predict_func <- function(obj, newx)
    {
      colnames(newx) <- paste0("X", 1:ncol(newx)) 
      predict(obj, newx) 
    }
  }

  resids <- residuals(fit_mean)
  fit_sigma <- ahead::mlf(resids^2, 
                          lags = 2L, 
                          fit_func=fit_func,
                          predict_func=predict_func,
                          type_pi=type_pi,
                          h=h, 
                          B = B) 
  z <- resids/sqrt(fitted(fit_sigma))
  fit_z <- ahead::conformalize(FUN=model_residuals, 
                               method=type_sim_conformalize,
                               nsim=B,
                               y=z, h=h)
  f <- as.numeric(fit_mean$mean) + matrix(fit_z$sims, ncol=B)*sqrt(pmax(matrix(fit_sigma$sims, ncol=B), 0)) # think about this clipping
  #f <- as.numeric(fit_mean$mean) + matrix(fit_z$sims, ncol=B)*sqrt(matrix(fit_sigma$sims, ncol=B))
  mean_f <- rowMeans(f)
  alpha <- 1 - level/100
  lower_bound <- apply(f, 1, function(x) quantile(x, probs=alpha/2))
  upper_bound <- apply(f, 1, function(x) quantile(x, probs=1-alpha/2))
  out <- list()
  out$x <- x
  out$level <- level
  out$resids_t_test <- t.test(resids, conf.level = level/100)
  out$resids_kpss_test <- suppressWarnings(tseries::kpss.test(resids))
  out$sims <- ts(f, start = start_preds, 
                 frequency = freq_x)
  if (!is.null(ml_method))
  {
    out$method <- paste0(ml_method, "ARCH")
  } else {
    out$method <- "ML-ARCH"
  }  
  out$model <- list(fit_mean, fit_sigma, fit_z)
  out$mean <- ts(mean_f, 
                 start = start_preds, 
                 frequency = freq_x)
  out$lower <- ts(lower_bound, 
                  start = start_preds, 
                  frequency = freq_x)
  out$upper <- ts(upper_bound, 
                  start = start_preds, 
                  frequency = freq_x)
  return(structure(out, class = "forecast"))
}