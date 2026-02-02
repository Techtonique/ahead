#' Ridge Regression Forecasting with Attention-based Context Vectors
#'
#' This function performs ridge regression forecasting using attention-based
#' context vectors as external regressors. Context vectors are computed from
#' the training data using various attention mechanisms and then used to
#' enhance the forecast model.
#'
#' @param y A multivariate time series object.
#' @param h Integer. Forecast horizon. Defaults to 5.
#' @param attention_type String specifying the type of attention mechanism.
#'   Options are: "cosine", "exponential", "dot_product", "scaled_dot_product",
#'   "gaussian", "linear", "value_based", "hybrid", "parametric". 
#'   Default is "exponential".
#' @param window_size Integer parameter for window size (applicable for 
#'   "cosine" attention). Defaults to 3.
#' @param decay_factor Double for decay factor (applicable for "exponential" 
#'   and "hybrid" attention). Defaults to 5.0.
#' @param temperature Double for temperature (applicable for 
#'   "scaled_dot_product" attention). Defaults to 1.0.
#' @param sigma Double for sigma (applicable for "gaussian" attention). 
#'   Defaults to 1.0.
#' @param sensitivity Double for sensitivity (applicable for "value_based" or 
#'   "hybrid" attention). Defaults to 1.0.
#' @param alpha Double for alpha (applicable for "parametric" attention). 
#'   Defaults to 0.5.
#' @param beta Double for beta (applicable for "parametric" attention). 
#'   Defaults to 0.5.
#' @param ... Additional arguments passed to \code{\link{ahead::ridge2f}}  
#'   (e.g., lags, lambda_1, lambda_2, nb_hidden, etc.).
#'
#' @return An object returned by \code{ahead::ridge2f} for the forecast,
#'   typically a list including \code{mean} and prediction intervals.
#'
#' @details
#' This approach allows the model to leverage temporal dependencies captured
#' by attention mechanisms, potentially improving forecast accuracy by
#' incorporating weighted historical information.
#'
#' @examples
#' 
#' plot(contextridge2f(AirPassengers, h = 15, lags = 15, attention_type = "exponential"))
#' 
#' plot(contextridge2f(fdeaths, h = 20, lags = 15, attention_type = "exponential"))
#' 
#' @seealso \code{\link{ahead::computeattention}}
#' 
#' @export
contextridge2f <- function(y,
                           h = 5L,
                           attention_type = "exponential",
                           window_size = 3,
                           decay_factor = 5.0,
                           temperature = 1.0,
                           sigma = 1.0,
                           sensitivity = 1.0,
                           alpha = 0.5,
                           beta = 0.5,
                           ...)
{
  if (is.null(dim(y)))
  {
    ctx_result <- ahead::computeattention(
      series = y,
      attention_type = attention_type,
      window_size = window_size,
      decay_factor = decay_factor,
      temperature = temperature,
      sigma = sigma,
      sensitivity = sensitivity,
      alpha = alpha,
      beta = beta
    ) 
    
    res<- ahead::ridge2f(
      y = y,
      h = h,
      xreg = ctx_result$context_vectors,
      ...
    )
    
    res$context_vectors <- ctx_result$context_vectors
    
  } else {
    ctx_result <- sapply(1:ncol(y), function(i) ahead::computeattention(
      series = y[,i],
      attention_type = attention_type,
      window_size = window_size,
      decay_factor = decay_factor,
      temperature = temperature,
      sigma = sigma,
      sensitivity = sensitivity,
      alpha = alpha,
      beta = beta
    )$context_vectors)
    res<- ahead::ridge2f(
      y = y,
      h = h,
      xreg = ctx_result,
      ...
    )
    res$context_vectors <- ctx_result
  }
  
  return(res)
}