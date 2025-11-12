#' @title Generalized Linear Model Theta Forecast and not only
#' @description This function implements the Theta method using a Generalized Linear Model (GLM)
#' @param y The time series data
#' @param h The number of periods to forecast
#' @param level The confidence level for the forecast intervals
#' @param fit_func The function to use for fitting the GLM (and other models)
#' @param predict_func The function to use for predict the other models (not GLMs)
#' @param fan Logical flag for fan plot
#' @param x The time series data
#' @param attention Logical flag for using attention mechanism
#' @param attention_type The type of attention mechanism to use
#' @param attention_method The method for computing attention-based adjustments: "heuristic" (original) or "historical" (new evidence-based)
#' @param scale_ctxt Scaling coefficient for context vector
#' @param historical_lookback Number of periods to consider for historical matching
#' @param max_adjustment Maximum allowed adjustment percentage (0-1)
#' @param B Number of bootstrap replications or number of simulations (alias: nsim)
#' @param nsim Alias for B
#' @param ... Additional arguments to pass to the fit_func
#' @return A forecast object
#' @export
glmthetaf <- function (
    y,
    h = ifelse(frequency(y) > 1, 2 * frequency(y), 10),
    level = 95L,
    fit_func = stats::glm,
    predict_func = predict, 
    fan = FALSE,
    x = y,
    type_pi = c("gaussian",
      "conformal-split", "conformal-surrogate",
      "conformal-bootstrap", "conformal-block-bootstrap",
      "conformal-kde", "conformal-fitdistr"),
    attention = TRUE,
    attention_type = c(
      "dot_product", "scaled_dot_product", "cosine", "exponential",
      "gaussian", "linear", "value_based", "hybrid", "parametric"
    ),
    attention_method = c("heuristic", "historical"),
    scale_ctxt = 1,
    historical_lookback = 50L,
    max_adjustment = 0.5,
    B = 250L,
    nsim = B,
    ...
)
{
  type_pi <- match.arg(type_pi)
  attention_type <- match.arg(attention_type)
  attention_method <- match.arg(attention_method)
  stopifnot(scale_ctxt > 0 && scale_ctxt <= 1)
  stopifnot(max_adjustment > 0 && max_adjustment <= 1)
  
  # --- Handle conformal variants via recursion ---
  if (grepl("conformal", type_pi)) {
    stopifnot(length(level) == 1)
    freq_x <- frequency(y)
    start_x <- start(y)
    start_preds <- tsp(y)[2] + 1 / freq_x
    # Split the training data
    y_train_calibration <- misc::splitts(y, split_prob = 0.5)
    y_train <- y_train_calibration$training
    y_calibration <- y_train_calibration$testing
    h_calibration <- length(y_calibration)
    # Predictions on calibration set
    y_pred_calibration <- ahead::glmthetaf(
      y_train,
      h = h_calibration,
      fit_func = fit_func,
      attention = attention,
      attention_type = attention_type,
      attention_method = attention_method,  
      type_pi = "gaussian",
      ...
    )$mean
    # Final fit and forecast on full calibration set
    fit_obj_train <- ahead::glmthetaf(
      y_calibration,
      h = h,
      fit_func = fit_func,
      attention = attention,
      attention_type = attention_type,
      attention_method = attention_method,
      type_pi = "gaussian",
      ...
    )
    
    fit_obj_train$method <- paste0(
      "Conformal (", gsub("conformal-", "", type_pi), ") Theta"
    )
    
    preds <- fit_obj_train$mean
    calibrated_residuals <- y_calibration - y_pred_calibration
    scaled_calib_resids <- scale(calibrated_residuals)
    sd_calibrated_residuals <- sd(calibrated_residuals)
    method <- gsub("conformal-", "", type_pi)
    
    if (method == "fitdistr") {
      simulate_function <- misc::fit_param_dist(as.numeric(scaled_calib_resids), verbose = FALSE)
      sims <- matrix(simulate_function(nsim * h), ncol = nsim, nrow = h)
      preds <- as.numeric(fit_obj_train$mean) + sims * sd_calibrated_residuals
      
      res <- list()
      res$level <- level
      res$x <- y_calibration
      res$method <- fit_obj_train$method
      start_preds <- start(fit_obj_train$mean)
      res$mean <- ts(rowMeans(preds), start = start_preds, frequency = freq_x)
      res$upper <- ts(apply(preds, 1, function(x)
        quantile(x, probs = 1 - (1 - level / 100) / 2)),
        start = start_preds, frequency = freq_x)
      res$lower <- ts(apply(preds, 1, function(x)
        quantile(x, probs = (1 - level / 100) / 2)),
        start = start_preds, frequency = freq_x)
      res$sims <- ts(preds, start = start_preds, frequency = freq_x)
      class(res) <- "forecast"
      return(res)
    }
    
    
    if (method == "block-bootstrap")
    {
      start_preds <- start(fit_obj_train$mean)
      sims <- ts(tseries::tsbootstrap(scaled_calib_resids, nb=nsim)[seq_len(h), ],
                 start = start_preds, 
                 frequency = freq_x) # nrow = h, ncol = nsim
      preds <- as.numeric(fit_obj_train$mean) + sims*sd_calibrated_residuals
      
      res <- list()
      res$level <- level 
      res$x <- y_calibration
      res$method <- fit_obj_train$method
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
      start_preds <- start(fit_obj_train$mean)
      sims <- ts(matrix(direct_sampling(scaled_calib_resids, 
                                        n = nsim*h, method=method), 
                        nrow = h, ncol = nsim), start = start_preds, 
                 frequency = freq_x)    
      preds <- as.numeric(fit_obj_train$mean) + sims*sd_calibrated_residuals
      
      res <- list()
      res$level <- level 
      res$x <- y_calibration
      res$method <- fit_obj_train$method
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
    
    if (method == "split") {
      quantile_absolute_residuals_conformal <- quantile(abs(scaled_calib_resids), 
                                                        probs = level/100)
      # Create output with proper time series attributes
      out <- list(
        mean = ts(preds, start = start_preds, frequency = freq_x),
        lower = ts(preds - sd_calibrated_residuals*quantile_absolute_residuals_conformal, 
                   start = start_preds, frequency = freq_x),
        upper = ts(preds + sd_calibrated_residuals*quantile_absolute_residuals_conformal, 
                   start = start_preds, frequency = freq_x),
        sims = NULL,
        x = y,
        level = level,
        method = fit_obj_train$method,
        residuals = ts(calibrated_residuals, 
                       start = start(y), 
                       frequency = freq_x)
      )
      return(structure(out, class = "forecast"))
    }
    
  } # end conformal section
  
  # --- Non-conformal: standard or attention-based Theta ---
  
  if (fan) {
    level <- seq(51, 99, by = 3)
  } else {
    if (min(level) > 0 && max(level) < 1) level <- 100 * level
    else if (min(level) < 0 || max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  
  n <- length(x)
  x <- as.ts(x)
  m <- frequency(x)
  seasonal <- FALSE
  if (m > 1 && !is.constant(x) && n > 2 * m) {
    r <- as.numeric(acf(x, lag.max = m, plot = FALSE)$acf)[-1]
    stat <- sqrt((1 + 2 * sum(r[-m]^2)) / n)
    seasonal <- (abs(r[m]) / stat > qnorm(0.95))
  }
  origx <- x
  
  if (seasonal) {
    decomp <- decompose(x, type = "multiplicative")
    if (any(abs(seasonal(decomp)) < 1e-04))
      warning("Seasonal indexes close to zero. Using non-seasonal Theta method")
    else
      x <- seasadj(decomp)
  }
  
  fcast <- ses(x, h = h)
  alpha <- pmax(1e-10, fcast$model$par["alpha"])
  time_idx <- 0:(n - 1)
  df <- data.frame(y = x, t = time_idx)
  
  # --- Case 1: No attention ---
  if (!attention) {
    tmp2 <- try(fit_func(y ~ t, data = df, ...), silent = TRUE)
    if (inherits(tmp2, "try-error")) {
      X <- matrix(time_idx, ncol = 1)
      tmp2 <- try(fit_func(x = X, y = as.numeric(x), ...), silent = TRUE)
      if (inherits(tmp2, "try-error"))
      {
        tmp2 <- try(fit_func(X, as.numeric(x), ...), silent = TRUE)
        if (inherits(tmp2, "try-error")) {
          stop("Unable to fit linear trend")
        }
      }
    }
    slope <- try(coef(tmp2)[2], silent = TRUE)
    if (inherits(slope, "try-error") | is.null(slope))
    {
      slope <- try(tmp2$coefficients[2], silent = TRUE)
      if (inherits(slope, "try-error") | is.null(slope))
      {
        slope <- estimate_theta_slope(fit_func = fit_func,
                                      predict_func = predict_func,
                                      y = y, ...)
      }
    }
    tmp2 <- slope / 2
    #misc::debug_print(slope)
    #misc::debug_print(tmp2)
    #misc::debug_print(fcast)
    fcast$mean <- fcast$mean + tmp2 * (0:(h - 1) + (1 - (1 - alpha)^n) / alpha)
  } else {
    tmp2 <- try(fit_func(y ~ t, data = df, ...), silent = TRUE)
    if (inherits(tmp2, "try-error")) {
      X <- matrix(time_idx, ncol = 1)
      tmp2 <- try(fit_func(x = X, y = x, ...), silent = TRUE)
      if (inherits(tmp2, "try-error"))
      {
        tmp2 <- try(fit_func(X, x, ...), silent = TRUE)
        if (inherits(tmp2, "try-error"))
        {
          stop("Unable to fit linear trend")
        }
      }
    }
    slope <- try(coef(tmp2)[2], silent = TRUE)
    if (inherits(slope, "try-error") | is.null(slope))
    {
      slope <- try(tmp2$coefficients[2], silent = TRUE)
      if (inherits(slope, "try-error") | is.null(slope))
      {
        slope <- estimate_theta_slope(fit_func = fit_func,
                                      predict_func = predict_func,
                                      y = y, ...)
      }
    }
    tmp2 <- slope / 2  # θ drift
    # Compute attention inside branch
    context_vectors <- ahead::computeattention(
      x, attention_type = attention_type, ...
    )$context_vectors
    last_context <- tail(context_vectors, 1)
    
    # --- NEW: Choose adjustment method ---
    if (attention_method == "heuristic") {
      # Original heuristic approach (maintains backward compatibility)
      context_adjusted_drift <- tmp2 * (1 + scale_ctxt * sign(last_context) *
                                          abs(last_context / pmax(mean(abs(y)), 1e-10) ))
    } else {
      # New historical evidence-based approach
      context_adjusted_drift <- historical_adjustment(
        base_drift = tmp2,
        current_context = last_context,
        full_series = y,
        historical_lookback = historical_lookback,
        max_adjustment = max_adjustment,
        scale_ctxt = scale_ctxt
      )
    }
    
    fcast$mean[1] <- fcast$mean[1] +
      context_adjusted_drift * ((1 - (1 - alpha)^n) / alpha)
    
    newx <- c(y, fcast$mean[1])
    for (i in 2:h) {
      context_vectors <- ahead::computeattention(
        newx, attention_type = attention_type, ...
      )$context_vectors
      last_context <- tail(context_vectors, 1)
      
      # Apply same method choice for recursive forecasts
      if (attention_method == "heuristic") {
        context_adjusted_drift <- tmp2 * (1 + scale_ctxt * sign(last_context) *
                                            abs(last_context / pmax(mean(abs(newx)), 1e-10)))
      } else {
        context_adjusted_drift <- historical_adjustment(
          base_drift = tmp2,
          current_context = last_context,
          full_series = newx,
          historical_lookback = min(historical_lookback, length(newx)),
          max_adjustment = max_adjustment,
          scale_ctxt = scale_ctxt
        )
      }
      
      fcast$mean[i] <- fcast$mean[i] +
        context_adjusted_drift * ((i - 1) + (1 - (1 - alpha)^n) / alpha)
      newx <- c(newx, fcast$mean[i])
    }
  }
  
  # --- Post-processing ---
  if (seasonal) {
    fcast$mean <- fcast$mean * rep(tail(decomp$seasonal, m),
                                   trunc(1 + h / m))[1:h]
    fcast$fitted <- fcast$fitted * decomp$seasonal
  }
  fcast$residuals <- origx - fcast$fitted
  fcast.se <- sqrt(fcast$model$sigma2) * sqrt((0:(h - 1)) * alpha^2 + 1)
  nconf <- length(level)
  fcast$lower <- fcast$upper <- ts(matrix(NA, nrow = h, ncol = nconf))
  tsp(fcast$lower) <- tsp(fcast$upper) <- tsp(fcast$mean)
  for (i in 1:nconf) {
    zt <- -qnorm(0.5 - level[i] / 200)
    fcast$lower[, i] <- fcast$mean - zt * fcast.se
    fcast$upper[, i] <- fcast$mean + zt * fcast.se
  }
  
  fcast$x <- origx
  fcast$level <- level
  fcast$method <- paste0("ML+Theta", if(attention) paste0(" with ", attention_method, " attention") else "")
  fcast$model <- list(alpha = alpha, drift = tmp2, sigma = fcast$model$sigma2,
                      attention_method = if(attention) attention_method else NULL)
  fcast$model$call <- match.call()
  return(fcast)
}

historical_adjustment <- function(base_drift, current_context, full_series, 
                                  historical_lookback = 50, max_adjustment = 0.5,
                                  scale_ctxt = 1) {
  n <- length(full_series)
  # Ensure we have enough history for meaningful comparison
  if (n < 10) {
    warning("Insufficient data for historical adjustment. Using base drift.")
    return(base_drift)
  }
  # Use available history, up to specified lookback
  lookback_window <- max(10, min(historical_lookback, n - 5))
  historical_data <- tail(full_series, lookback_window)
  # Extract context features from current window
  current_features <- c(
    level = mean(historical_data),
    volatility = sd(historical_data),
    trend = if(length(historical_data) > 1) {
      coef(lm(historical_data ~ seq_along(historical_data)))[2]
    } else 0,
    recent_momentum = if(length(historical_data) > 1) {
      tail(historical_data, 1) - head(historical_data, 1)
    } else 0
  )
  # Find similar historical patterns using sliding windows
  window_size <- length(historical_data)
  similarities <- numeric()
  optimal_adjustments <- numeric()
  
  if (n > 2 * window_size) {
    for (i in 1:(n - 2 * window_size + 1)) {
      # Extract historical window and subsequent window
      hist_window <- full_series[i:(i + window_size - 1)]
      subsequent_window <- full_series[(i + window_size):min(i + 2 * window_size - 1, n)]
      # Compute features for historical window
      hist_features <- c(
        level = mean(hist_window),
        volatility = sd(hist_window),
        trend = coef(lm(hist_window ~ seq_along(hist_window)))[2],
        recent_momentum = tail(hist_window, 1) - head(hist_window, 1)
      )
      # Compute similarity (Euclidean distance on normalized features)
      similarity <- 1 / (1 + sqrt(sum((current_features - hist_features)^2)))
      similarities <- c(similarities, similarity)
      # Compute what adjustment would have been optimal
      actual_subsequent_trend <- if(length(subsequent_window) > 1) {
        coef(lm(subsequent_window ~ seq_along(subsequent_window)))[2]
      } else 0
      # Optimal adjustment = (actual_trend - base_drift) / base_drift
      optimal_adj <- (actual_subsequent_trend - base_drift) / abs(base_drift + 1e-10)
      optimal_adjustments <- c(optimal_adjustments, optimal_adj)
    }
  }
  # If we found historical matches, compute weighted adjustment
  if (length(similarities) > 0 && max(similarities) > 0.1) {
    # Use top-k most similar historical patterns
    k <- min(5, length(similarities))
    top_indices <- order(similarities, decreasing = TRUE)[1:k]
    top_similarities <- similarities[top_indices]
    top_adjustments <- optimal_adjustments[top_indices]
    # Weight adjustments by similarity
    weights <- top_similarities / sum(top_similarities)
    suggested_adjustment <- sum(weights * top_adjustments)
    # Confidence based on similarity consistency
    confidence <- 1 - (sd(top_adjustments) / (abs(mean(top_adjustments)) + 1e-10))
    confidence <- max(0, min(1, confidence))
    # Apply scaling and bounds
    final_adjustment <- scale_ctxt * suggested_adjustment * confidence
    bounded_adjustment <- sign(final_adjustment) * 
      min(abs(final_adjustment), max_adjustment)
    adjusted_drift <- base_drift * (1 + bounded_adjustment)
  } else {
    # No meaningful historical matches found
    warning("No similar historical patterns found. Using base drift with context scaling.")
    # Fallback to scaled version of original heuristic
    context_effect <- scale_ctxt * sign(current_context) * 
      min(abs(current_context / (sd(full_series) + 1e-10)), 1)
    adjusted_drift <- base_drift * (1 + context_effect * max_adjustment)
  }
  
  return(adjusted_drift)
}