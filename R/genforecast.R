  #' 
#' @title Generic Forecasting Function (Unified interface)
#' 
#' @description
#' 
#' This function allows to call any function "of class \code{forecast}" in a unified way.
#' 
#' @param FUN A forecasting function.
#' @param y A time series (\code{ts} object or vector).
#' @param h Forecasting horizon.
#' @param level The confidence level.
#' @param ... Additional arguments to be passed to the forecasting function.
#' 
#' @return An object of class \code{forecast}.
#' 
#' @examples
#' 
#' y <- fdeaths 
#' h <- 25L
#' plot(ahead::genericforecast(FUN=forecast::thetaf, y, h))
#' plot(ahead::genericforecast(FUN=ahead::dynrmf, y, h))
#' plot(ahead::genericforecast(FUN=forecast::tbats, y=y, h=h, use.box.cox = TRUE, use.trend=FALSE))
#' plot(ahead::genericforecast(FUN=forecast::auto.arima, y, h)) 
#'
#' @export
#' 
genericforecast <- function(FUN, y, h, level=95, ...)
{
  obj <- try(do.call(what=FUN, args=list(y = y, h=h, level=level, ...)), 
             silent = TRUE) # forecast::thetaf e.g 

  if (inherits(obj, "forecast"))
  {
    return(obj)
  }

  if (inherits(obj, "try-error"))
  {
    obj <- try(do.call(what=FUN, args=list(y = y, ...)), 
               silent = TRUE) # forecast::auto.arima or forecast::ets e.g 
    if (inherits(obj, "try-error"))
    {
      obj <- try(do.call(what=FUN, args=list(x = y, ...)), 
               silent = TRUE) # Holtwinters e.g
    }
  }

  res <- try(forecast::forecast(obj, h=h, level=level, ...), 
             silent = TRUE)

  if (inherits(res, "try-error"))
  {
    res <- do.call(what = predict, 
                   args = list(object=obj, n.ahead=h, prediction.interval = TRUE, level = level, ...))
  }
  return(res)
}