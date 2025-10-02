#' Transformer Forecasting for Financial Returns
#'
#' @param y Univariate time series of returns
#' @param h Forecast horizon
#' @param level Confidence levels for prediction intervals
#' @param sequence_size Lookback period for transformer
#' @param ... Additional parameters to TRANSFORMER function
#' 
#' @export 
#' 
#' @return Object of class 'forecast'
#' 
transformerreturnsf <- function(y, h = 10, level = c(80, 95), 
                                sequence_size = 20, num_heads = 4,
                                head_size = 32, num_transformer_blocks = 2,
                                epochs = 50, ...) {
  if (!requireNamespace("transformerForecasting", quietly = TRUE)) {
    stop("transformerForecasting package required")
  }
  # Input validation for returns data
  if (is.ts(y)) {
    freq <- frequency(y)
    start <- tsp(y)[1]
  } else {
    freq <- 1
    start <- 1
  }
  # Prepare data frame for transformerForecasting
  n <- length(y)
  df <- data.frame(
    Date = 1:n,
    return = as.numeric(y)
  )
  # Fit transformer model to historical data
  transformer_result <- suppressWarnings(transformerForecasting::TRANSFORMER(
    df, 
    study_variable = "return",
    sequence_size = sequence_size,
    num_heads = num_heads,
    head_size = head_size,
    num_transformer_blocks = num_transformer_blocks,
    epochs = epochs,
    ...
  ))
  # Extract predictions (these are one-step-ahead forecasts)
  predictions <- transformer_result$PREDICTIONS
  # For future forecasts, use recursive strategy with the last sequence
  future_forecasts <- numeric(h)
  # Start with the full historical data
  current_data <- as.numeric(y)
  
  for(i in 1:h) {
    # Create temporary data frame with current available data
    temp_n <- length(current_data)
    temp_df <- data.frame(
      Date = 1:temp_n,
      return = current_data
    )
    
    # Get transformer prediction
    temp_result <- suppressWarnings(transformerForecasting::TRANSFORMER(
      temp_df,
      study_variable = "return",
      sequence_size = min(sequence_size, temp_n),  # Ensure sequence_size <= data length
      num_heads = num_heads,
      head_size = head_size,
      num_transformer_blocks = num_transformer_blocks,
      epochs = 1#,  # Single epoch for prediction
      #verbose = 0
    ))
    # Take the last prediction as the next forecast
    next_forecast <- tail(temp_result$PREDICTIONS, 1)
    future_forecasts[i] <- next_forecast
    # Append the forecast to current data for next iteration
    current_data <- c(current_data, next_forecast)
  }
  # Calculate fitted values (remove NA from beginning due to sequence_size)
  fitted_vals <- c(rep(NA, sequence_size), 
                   predictions[(sequence_size + 1):n])
  # Calculate residuals
  residuals_vals <- as.numeric(y) - fitted_vals
  # Generate prediction intervals using residual distribution
  valid_residuals <- residuals_vals[!is.na(residuals_vals)]
  residual_sd <- sd(valid_residuals)
  # Create prediction intervals
  upper <- lower <- list()
  for (lvl in level) {
    z <- qnorm(1 - (1 - lvl/100)/2)
    # Scale error margin with sqrt(horizon) for increasing uncertainty
    error_margin <- z * residual_sd * sqrt(seq_len(h))
    upper[[as.character(lvl)]] <- future_forecasts + error_margin
    lower[[as.character(lvl)]] <- future_forecasts - error_margin
  }
  # Return forecast object
  result <- list(
    method = "Transformer Returns",
    model = transformer_result,
    mean = ts(future_forecasts, frequency = freq, 
              start = c(tsp(y)[2] + 1/freq, 1)),
    x = y,
    fitted = ts(fitted_vals, frequency = freq, start = start),
    residuals = ts(residuals_vals, frequency = freq, start = start),
    level = level,
    upper = do.call(cbind, upper),
    lower = do.call(cbind, lower),
    sequence_size = sequence_size
  )
  class(result) <- "forecast"
  return(result)
}

# Helper function for NULL coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b