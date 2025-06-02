#' @export 
mlarchf <- function(y, h=10L,
                    mean_model=forecast::auto.arima,
                    model_residuals=forecast::thetaf,
                    fit_func=ahead::ridge,
                    predict_func=predict,
                    level=95,
                    B=250L)
{
  x <- y
  last_price <- as.numeric(tail(y, 1))
  freq_x <- frequency(y)
  start_preds <- tsp(y)[2] + 1/freq_x
  fit_init_model <- ahead::genericforecast(mean_model, 
                                           y=y, h=h)
  resids_init_model <- residuals(fit_init_model)
  fit_sigma <- ahead::mlf(resids_init_model^2, 
                            lags = 2L, 
                            fit_func=fit_func,
                            predict_func=predict_func,
                            h=h, 
                            B = B) 
  z <- resids_init_model/sqrt(fitted(fit_sigma))
  fit_z <- ahead::conformalize(FUN=model_residuals, 
                               nsim=B,
                               y=z, h=h)
  f <- as.numeric(fit_init_model$mean) + matrix(fit_z$sims, ncol=B)*sqrt(pmax(matrix(fit_sigma$sims, ncol=B), 0))
  mean_f <- rowMeans(f)
  alpha <- 1 - level/100
  lower_bound <- apply(f, 1, function(x) quantile(x, probs=alpha/2))
  upper_bound <- apply(f, 1, function(x) quantile(x, probs=1-alpha/2))
  out <- list()
  out$x <- x
  out$level <- level
  out$t_test <- t.test(resids_init_model, conf.level = level/100)
  out$kpss_test <- tseries::kpss.test(resids_init_model)
  out$method <- paste0(method, "ARCH")
  out$model <- list(fit_init_model, fit_sigma, fit_z)
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