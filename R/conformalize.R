#' Conformalize a forecasting function 
#' 
#' @description 
#' 
#' This function allows to conformalize any forecasting function.
#' 
#' @param FUN A forecasting function.
#' @param y A time series (\code{ts} object or vector).
#' @param h Forecasting horizon.
#' @param level Confidence level.
#' @param method Method to be used for conformalization (simulation of calibrated residuals).
#' @param nsim Number of simulations.
#' @param block_size Block size for block-bootstrap.
#' @param seed Seed for reproducibility.
#' @param ... Additional arguments to be passed to the forecasting function.
#' 
#' @return An object of class \code{forecast}.
#' 
#' @examples
#' 
#' y <- fdeaths
#' h <- 25L
#' obj <- conformalize(FUN=forecast::ets, y, h); plot(obj)
#' obj <- conformalize(FUN=HoltWinters, y=y, h=h, seasonal = "mult"); plot(obj)
#' 
conformalize <- function(FUN, y, h, level=95,
                         method = c("block-bootstrap", "surrogate", "kde", "bootstrap"),
                         nsim = 100L, 
                         block_size = 5,
                         seed = 123L, 
                         ...)
{
  method <- match.arg(method)
  n <- length(y)
  freq_x <- frequency(y)
  half_n <- ceiling(n/2)
  idx_train <- seq_len(half_n)
  idx_calib <- setdiff(seq_len(n), idx_train)
  splitted_y <- misc::splitts(y)
  y_train <- splitted_y$training
  y_calib <- splitted_y$testing
  n_calib <- length(idx_calib)
  
  calib_resids <- ahead::genericforecast(y=y_train, h=n_calib, 
                              level=level, FUN=FUN, 
                              ...)$mean - y_calib
  scaled_calib_resids <- base::scale(calib_resids)
  xm <- attr(scaled_calib_resids, "scaled:center")
  xsd <- attr(scaled_calib_resids, "scaled:scale")
  obj_fcast <- ahead::genericforecast(y=y_calib, h=h, 
                          level=level, FUN=FUN, 
                          ...) # train on calibration set 

   if (method == "fitdistr")
  {                            
    simulate_function <- misc::fit_param_dist(as.numeric(scaled_calib_resids), 
                                              verbose = FALSE)
    sims <- matrix(simulate_function(nsim*h), ncol=nsim, nrow=h)
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- y_calib
    res$method <- paste0("conformalized ", obj_fcast$method)
    start_preds <- start(obj_fcast$mean)
    res$mean <- ts(rowMeans(preds), 
                   start = start_preds, 
                   frequency = freq_x)
    res$upper <- ts(apply(preds, 1, function(x)
      stats::quantile(x, probs = 1 - (1 - level / 100) / 2)), 
      start = start_preds, 
      frequency = freq_x)
    res$lower <- ts(apply(preds, 1, function(x)
      stats::quantile(x, probs = (1 - level / 100) / 2)), 
      start = start_preds, 
      frequency = freq_x)
    res$sims <- ts(preds, start = start_preds, 
                   frequency = freq_x)
    class(res) <- "forecast"
    return(res)
  }

  if (method == "block-bootstrap")
  {
    start_preds <- start(obj_fcast$mean)
    sims <- ts(sapply(1:nsim, function(i) mbb(matrix(scaled_calib_resids, 
                                                 ncol = 1), 
                                          n=h, 
                                          b=block_size, 
                                          seed=i+seed*100)),
             start = start_preds, 
             frequency = freq_x)
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- y_calib
    res$method <- paste0("conformalized ", obj_fcast$method)    
    res$mean <- ts(rowMeans(preds), 
                   start = start_preds, 
                   frequency = freq_x)
    res$upper <- ts(apply(preds, 1, function(x)
      stats::quantile(x, probs = 1 - (1 - level / 100) / 2)), 
      start = start_preds, 
      frequency = freq_x)
    res$lower <- ts(apply(preds, 1, function(x)
      stats::quantile(x, probs = (1 - level / 100) / 2)), 
      start = start_preds, 
      frequency = freq_x)
    res$sims <- ts(preds, start = start_preds, 
                   frequency = freq_x)
    class(res) <- "forecast"
    return(res)
  }

  if (method %in% c("kde", "surrogate", "bootstrap"))
  {
    start_preds <- start(obj_fcast$mean)
    sims <- ts(matrix(direct_sampling(scaled_calib_resids, 
    n = nsim*h, method=method), 
    nrow = h, ncol = nsim), start = start_preds, 
             frequency = freq_x)    
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- y_calib
    res$method <- paste0("conformalized ", obj_fcast$method)    
    res$mean <- ts(rowMeans(preds), 
                   start = start_preds, 
                   frequency = freq_x)
    res$upper <- ts(apply(preds, 1, function(x)
      stats::quantile(x, probs = 1 - (1 - level / 100) / 2)), 
      start = start_preds, 
      frequency = freq_x)
    res$lower <- ts(apply(preds, 1, function(x)
      stats::quantile(x, probs = (1 - level / 100) / 2)), 
      start = start_preds, 
      frequency = freq_x)
    res$sims <- ts(preds, start = start_preds, 
                   frequency = freq_x)
    class(res) <- "forecast"
    return(res)
  }
}