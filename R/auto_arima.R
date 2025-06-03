auto_arima <- function(y, h=5L, n_starts=10L, ...)
{
  penalty <- 1e6
  ndims <- 6
  objective <- function(xx)
  {
    model <- try(forecast::Arima(y=y, 
                             order = c(floor(xx[1]), 
                                       floor(xx[2]), 
                                       floor(xx[3])), 
                             seasonal = c(floor(xx[4]), 
                                          floor(xx[5]), 
                                          floor(xx[6])), 
                             method = "ML"), 
                 silent = TRUE)
    if (inherits(model, "try-error"))
    {
      return(.Machine$double.xmax)
    }
    ans <- try(model$aicc + penalty*(suppressWarnings(tseries::kpss.test(model$residuals))$p.value < 0.05),
               silent = TRUE)
    cond <- inherits(ans, "try-error")
    return(ans*(!cond)+.Machine$double.xmax*cond)
  }
  objective <- memoise::memoise(objective)
  
  res <- vector("list", length = n_starts)
  lower_bounds <- rep(0, ndims)
  upper_bounds <- rep(5, ndims)
  upper_bounds[2] <- 2
  upper_bounds[4] <- 2
  upper_bounds[6] <- 12
  pb <- utils::txtProgressBar(min=1, max=n_starts, style=3L)
  for (i in seq_len(n_starts))
  {
    set.seed(123 + i * 100)
    res[[i]] <- stats::nlminb(start = floor(3*runif(ndims)), 
                              objective = objective, 
                              lower = lower_bounds, 
                              upper = upper_bounds)
    if (i >= 2)
    {
      if (res[[i]]$objective <= res[[i-1]]$objective)
      {
        ans <- res[[i]]
      } 
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  return(forecast::Arima(y=y, 
                         order = c(floor(ans$par[1]), 
                                   floor(ans$par[2]), 
                                   floor(ans$par[3])), 
                         seasonal = c(floor(ans$par[4]), 
                                       floor(ans$par[5]), 
                                       floor(ans$par[6])), 
                         method = "ML"))
}

