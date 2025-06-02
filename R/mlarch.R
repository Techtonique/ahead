#' @export 
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
  x <- y
  last_price <- as.numeric(tail(y, 1))
  freq_x <- frequency(y)
  start_preds <- tsp(y)[2] + 1/freq_x
  type_pi <- match.arg(type_pi)
  type_sim_conformalize <- match.arg(type_sim_conformalize)
  fit_mean <- ahead::genericforecast(mean_model, 
                                           y=y, h=h)
  # conformalize mean model maybe?
  if (!is.null(method))
  {
      fit_func <- function(x, y, method = "ranger", ...)
    {
      df <- data.frame(y=y, as.matrix(x)) # naming of columns is mandatory for `predict`
      colnames(df) <- c("y", paste0("X", 1:ncol(x)))
      caret::train(y ~ ., data=df,
                  method = method,
                  trControl=caret::trainControl(method = "none"))
    }

    predict_func <- function(obj, newx)
    {
      colnames(newx) <- paste0("X", 1:ncol(newx)) # mandatory, linked to df in fit_func
      predict(obj, newx) # only accepts a named newx
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
  f <- as.numeric(fit_mean$mean) + matrix(fit_z$sims, ncol=B)*sqrt(pmax(matrix(fit_sigma$sims, ncol=B), 0))
  mean_f <- rowMeans(f)
  alpha <- 1 - level/100
  lower_bound <- apply(f, 1, function(x) quantile(x, probs=alpha/2))
  upper_bound <- apply(f, 1, function(x) quantile(x, probs=1-alpha/2))
  out <- list()
  out$x <- x
  out$level <- level
  out$t_test <- t.test(resids, conf.level = level/100)
  out$kpss_test <- tseries::kpss.test(resids)
  out$method <- paste0(method, "ARCH")
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