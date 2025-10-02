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
transformerreturnsf <- function(y, h = 10, level = 95, 
                                sequence_size = 20, num_heads = 4,
                                head_size = 32, num_transformer_blocks = 2,
                                epochs = 50, n_bootstrap = 250, 
                                show_progress = TRUE, seed=123, ...) {
  if (!requireNamespace("transformerForecasting", quietly = TRUE)) {
    stop("transformerForecasting package required. Please install it.")
  }
  set.seed(seed)
  # Input validation for returns data
  if (is.ts(y)) {
    freq <- frequency(y)
    start <- tsp(y)[1]
  } else {
    stop("Input y must be a time series (ts) object.")
  }
  # Function to apply the transformer forecasting recursively
  forecast_recursive <- function(y, h, sequence_size, num_heads, 
                                 head_size, num_transformer_blocks, 
                                 epochs, ...) {
    n <- length(y)
    df <- data.frame(
      Date = time(y),
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
    # Initialize the future forecast array
    future_forecasts <- numeric(h)
    current_data <- as.numeric(y)
    # For recursive forecasting
    for (i in 1:h) {
      temp_n <- length(current_data)
      temp_df <- data.frame(
        Date = 1:temp_n,
        return = current_data
      )
      temp_result <- suppressWarnings(transformerForecasting::TRANSFORMER(
        temp_df,
        study_variable = "return",
        sequence_size = min(sequence_size, temp_n),
        num_heads = num_heads,
        head_size = head_size,
        num_transformer_blocks = num_transformer_blocks,
        epochs = 1  # Single epoch for prediction
      ))
      # Get the next forecast
      next_forecast <- tail(temp_result$PREDICTIONS, 1)
      future_forecasts[i] <- next_forecast
      current_data <- c(current_data, next_forecast)
    }
    return(future_forecasts)
  }
  forecast_recursive <- compiler::cmpfun(forecast_recursive)
  
  # Bootstrap procedure
  boot_samples <- sapply(seq_len(n_bootstrap), function (i) ahead::rmultivariate(y, n = h))
  bootstrap_forecasts <- matrix(0, nrow=h, ncol=n_bootstrap)
  if (show_progress)
    pb <- utils::txtProgressBar(max=n_bootstrap, style = 3)
  for (i in seq_len(n_bootstrap)) 
  {
    # Call the forecasting function recursively on the bootstrapped sample
    bootstrap_forecasts[,i] <- forecast_recursive(boot_samples[,i], 
                                                  h, sequence_size, 
                                                  num_heads, head_size, 
                                                  num_transformer_blocks, epochs, ...)
    if (show_progress)
      utils::setTxtProgressBar(pb, i)
  }
  if (show_progress)
    close(pb)
  # Calculate prediction intervals (e.g., 5th and 95th percentiles)
  forecast_upper <- apply(bootstrap_forecasts, 1, function(x) quantile(x, level/100))
  forecast_lower <- apply(bootstrap_forecasts, 1, function(x) quantile(x, 1 - level/100))
  # Calculate the mean forecast (e.g., the average of the bootstrapped forecasts)
  forecast_mean <- apply(bootstrap_forecasts, 1, mean)
  # Return forecast object with updated intervals
  result <- list(
    method = "Transformer with Bootstrapping",
    model = NULL,  # No model stored here since we're bootstrapping
    mean = ts(forecast_mean, frequency = freq, start = c(tsp(y)[2] + 1 / freq, 1)),
    x = y,
    upper = ts(forecast_upper, frequency = freq, start = c(tsp(y)[2] + 1 / freq, 1)),
    lower = ts(forecast_lower, frequency = freq, start = c(tsp(y)[2] + 1 / freq, 1)),
    sims = ts(bootstrap_forecasts, frequency = freq, start = c(tsp(y)[2] + 1 / freq, 1)),
    level = level,
    sequence_size = sequence_size
  )
  class(result) <- "forecast"
  return(result)
}

# Helper function for NULL coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b