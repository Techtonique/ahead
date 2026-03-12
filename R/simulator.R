#' simulate from a forecasting function 
#' 
#' @description 
#' 
#' This function allows to obtain predictive simulations from 
#' any forecasting function.
#' 
#' @param FUN A forecasting function.
#' @param y A time series (\code{ts} object or vector).
#' @param h Forecasting horizon.
#' @param level Confidence level.
#' @param method Method to be used for the simulation of in-sample residuals.
#' @param nsim Number of simulations.
#' @param block_size Block size for block-bootstrap.
#' @param seed Seed for reproducibility.
#' @param B Alias for \code{nsim}
#' @param ... Additional arguments to be passed to the forecasting function.
#' 
#' @return An object of class \code{forecast}, with a slot \code{sims} 
#' 
#' @examples
#' 
#' y <- fdeaths
#' h <- 25L
#' obj <- simulator(FUN=forecast::thetaf, y, h); plot(obj)
#' obj <- simulator(FUN=forecast::auto.arima, y, h); plot(obj)
#' obj <- simulator(FUN=forecast::ets, y, h); plot(obj)
#'  
simulator <- function(FUN, y, h, level=95,
                      method = c("block-bootstrap", "surrogate", 
                                "kde", "bootstrap", "fitdistr", 
                                "meboot"),
                     nsim = 100L, 
                     block_size = 5,
                     seed = 123L, 
                     B = NULL,
                     ...)
{
  method <- match.arg(method)
  n <- length(y)
  freq_x <- frequency(y)
  if (!is.null(B))
    nsim <- B 
  obj_fcast <- ahead::genericforecast(y=y, h=h, 
                                level=level, FUN=FUN, 
                                ...)
  training_resids <- residuals(obj_fcast)
  scaled_training_resids <- base::scale(training_resids)
  xm <- attr(scaled_training_resids, "scaled:center")
  xsd <- attr(scaled_training_resids, "scaled:scale")

   if (method == "fitdistr")
  {                            
    simulate_function <- misc::fit_param_dist(as.numeric(scaled_training_resids), 
                                              verbose = FALSE)
    sims <- matrix(simulate_function(nsim*h), ncol=nsim, nrow=h)
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- 
    res$fitted <- fitted(obj_fcast)
    res$method <- paste0("simulated ", obj_fcast$method)
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
    sims <- ts(sapply(1:nsim, function(i) mbb(matrix(scaled_training_resids, 
                                                 ncol = 1), 
                                          n=h, 
                                          b=block_size, 
                                          seed=i+seed*100)),
             start = start_preds, 
             frequency = freq_x)
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- y
    res$fitted <- fitted(obj_fcast)
    res$method <- paste0("simulated ", obj_fcast$method)    
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
    sims <- ts(matrix(direct_sampling(scaled_training_resids, 
    n = nsim*h, method=method), 
    nrow = h, ncol = nsim), start = start_preds, 
             frequency = freq_x)    
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- y
    res$fitted <- fitted(obj_fcast)
    res$method <- paste0("simulated ", obj_fcast$method)    
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
  
  if (method == "meboot")
  {
    start_preds <- start(obj_fcast$mean)
    boots <- ahead::meboot(scaled_training_resids, 
                           reps = nsim)
    sims <- ts(boots$ensemble[seq_len(h), ], 
               start = start_preds, 
               frequency = freq_x)    
    preds <- as.numeric(obj_fcast$mean) + sims*xsd
    
    res <- list()
    res$level <- level 
    res$x <- y
    res$fitted <- fitted(obj_fcast)
    res$method <- paste0("simulated ", obj_fcast$method)    
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