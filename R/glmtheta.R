#' @title Generalized Linear Model Theta Forecast
#' @description This function implements the Theta method using a Generalized Linear Model (GLM)
#' @param y The time series data
#' @param h The number of periods to forecast
#' @param level The confidence level for the forecast intervals
#' @param fit_func The function to use for fitting the GLM
#' @param fan Logical flag for fan plot
#' @param x The time series data
#' @param attention Logical flag for using attention mechanism
#' @param scale_ctxt Scaling coefficient for context vector 
#' @param ... Additional arguments to pass to the fit_func
#' @return A forecast object
#' @export
glmthetaf <- function (y, 
                       h = ifelse(frequency(y) > 1, 2 * frequency(y), 10), 
                       level = 95, 
                       fit_func = stats::glm, 
                       fan = FALSE, 
                       x = y, 
                       type_pi = c(
                         "conformal-split",
                         "gaussian"
                       ),
                       attention = TRUE, 
                       scale_ctxt = 1,
                       method = c("adj", "base"),
                       ...) 
{
  method <- match.arg(method) 
  type_pi <- match.arg(type_pi)
  stopifnot(scale_ctxt > 0 && scale_ctxt <= 1)
  if (grepl("conformal", type_pi) >= 1) # not conformal
  {
    stopifnot(length(level) == 1)
    freq_x <- frequency(y)
    start_x <- start(y)
    start_preds <- tsp(y)[2] + 1 / freq_x
    # Split the training data
    y_train_calibration <- splitts(y, split_prob=0.5)
    y_train <- y_train_calibration$training 
    y_calibration <- y_train_calibration$testing
    h_calibration <- length(y_calibration)
    # Get predictions on calibration set
    y_pred_calibration <- ahead::glmthetaf(
      y_train,
      h = h_calibration,
      fit_func = fit_func, 
      attention = attention,
      method = method,
      type_pi = "gaussian"
    )$mean
    # Final fit and forecast on full calibration set
    fit_obj_train <- ahead::glmthetaf(
      y_calibration, 
      h = h,
      fit_func = fit_func,
      attention = attention,
      method = method,
      type_pi = "gaussian" 
    )
    preds <- fit_obj_train$mean
    calibrated_residuals <- y_calibration - y_pred_calibration
    scaled_calibrated_residuals <- base::scale(calibrated_residuals)
    sd_calibrated_residuals <- sd(calibrated_residuals)
  
    if (type_pi == "conformal-split") {
        quantile_absolute_residuals_conformal <- quantile(abs(scaled_calibrated_residuals), 
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
        method = "ConformalTheta",
        residuals = ts(calibrated_residuals, 
                       start = start(y), 
                       frequency = freq_x)
      )
      return(structure(out, class = "forecast"))
    }
  }
  
  if (fan) {
    level <- seq(51, 99, by = 3)
  }
  else {
    if (min(level) > 0 && max(level) < 1) {
      level <- 100 * level
    }
    else if (min(level) < 0 || max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  n <- length(x)
  x <- as.ts(x)
  m <- frequency(x)
  if (m > 1 && !is.constant(x) && n > 2 * m) {
    r <- as.numeric(acf(x, lag.max = m, plot = FALSE)$acf)[-1]
    stat <- sqrt((1 + 2 * sum(r[-m]^2))/n)
    seasonal <- (abs(r[m])/stat > qnorm(0.95))
  }
  else {
    seasonal <- FALSE
  }
  origx <- x
  if (seasonal) {
    decomp <- decompose(x, type = "multiplicative")
    if (any(abs(seasonal(decomp)) < 1e-04)) {
      warning("Seasonal indexes close to zero. Using non-seasonal Theta method")
    }
    else {
      x <- seasadj(decomp)
    }
  }
  
  fcast <- ses(x, h = h) 
  alpha <- pmax(1e-10, fcast$model$par["alpha"])
  
  # Try formula interface first, then matrix interface
  time_idx <- 0:(n - 1)
  df <- data.frame(y = x, t = time_idx)

  if (attention == FALSE) {
    tmp2 <- try(fit_func(y ~ t, data = df, ...), silent = TRUE)
    if (inherits(tmp2, "try-error")) {
      # For matrix interface, include intercept term for methods like glmnet
      X <- matrix(time_idx, nrow=1)
      tmp2 <- try(fit_func(x = X, y = as.numeric(x), ...), silent = TRUE)
      if (inherits(tmp2, "try-error")) {
        stop("Unable to fit linear trend")
      }
    }
    # Extract coefficient using generic method first, then fallback
    slope <- try(coef(tmp2)[2], silent = TRUE)
    if (inherits(slope, "try-error")) {
      slope <- try(tmp2$coefficients[2], silent = TRUE)
      if (inherits(slope, "try-error")) {
        slope <- try(tmp2$coef[2], silent = TRUE)
        if (inherits(slope, "try-error")) {
          stop("Unable to extract slope coefficient")
        }
      }
    }
    tmp2 <- slope/2  # Divide by 2 as per theta method
    fcast$mean <- fcast$mean + tmp2 * (0:(h - 1) + (1 - (1 - 
                                                           alpha)^n)/alpha)
    if (seasonal) {
      fcast$mean <- fcast$mean * rep(tail(decomp$seasonal, 
                                          m), trunc(1 + h/m))[1:h]
      fcast$fitted <- fcast$fitted * decomp$seasonal
    }
    fcast$residuals <- origx - fcast$fitted
    fcast.se <- sqrt(fcast$model$sigma2) * sqrt((0:(h - 1)) * 
                                                  alpha^2 + 1)
    nconf <- length(level)
    fcast$lower <- fcast$upper <- ts(matrix(NA, nrow = h, ncol = nconf))
    tsp(fcast$lower) <- tsp(fcast$upper) <- tsp(fcast$mean)
    for (i in 1:nconf) {
      zt <- -qnorm(0.5 - level[i]/200)
      fcast$lower[, i] <- fcast$mean - zt * fcast.se
      fcast$upper[, i] <- fcast$mean + zt * fcast.se
    }
    fcast$x <- origx
    fcast$level <- level
    fcast$method <- "Theta"
    fcast$model <- list(alpha = alpha, drift = tmp2, sigma = fcast$model$sigma2)
    fcast$model$call <- match.call()
    return(fcast)
  } else { # attention is TRUE
      
      if (method == "base")
      {
        #misc::debug_print(method)
        df <- cbind.data.frame(df, ctx = ahead::computeattention(x)$context_vectors)
        tmp2 <- try(fit_func(y ~ ., data = df, ...), silent = TRUE)
        if (inherits(tmp2, "try-error")) {
          X <- cbind.data.frame(time_idx, df$ctx)
          tmp2 <- try(fit_func(x = as.matrix(X), y = as.numeric(x), ...), silent = TRUE)
          if (inherits(tmp2, "try-error")) {
            stop("Unable to fit linear trend with attention")
          }
        }
      } else { # method == "adj": same as attention == FALSE to avoid adjusting the trend twice
        #misc::debug_print(method)
        tmp2 <- try(fit_func(y ~ t, data = df, ...), silent = TRUE)
        if (inherits(tmp2, "try-error")) {
          # For matrix interface, include intercept term for methods like glmnet
          X <- matrix(time_idx, nrow=1)
          tmp2 <- try(fit_func(x = X, y = as.numeric(x), ...), silent = TRUE)
          if (inherits(tmp2, "try-error")) {
            stop("Unable to fit linear trend")
          }
        }
        # Extract coefficient using generic method first, then fallback
        slope <- try(coef(tmp2)[2], silent = TRUE)
        if (inherits(slope, "try-error")) {
          slope <- try(tmp2$coefficients[2], silent = TRUE)
          if (inherits(slope, "try-error")) {
            slope <- try(tmp2$coef[2], silent = TRUE)
            if (inherits(slope, "try-error")) {
              stop("Unable to extract slope coefficient")
            }
          }
        }
        tmp2 <- slope/2  # Divide by 2 as per theta method
      }
    }
   
  if (method == 'base')
  {
    context_vectors <- log(ahead::computeattention(x)$context_vectors)
    last_context <- tail(context_vectors, 1)
    time_coef <- try(coef(tmp2)[2], silent = TRUE)
    ctx_coef <- try(coef(tmp2)[3], silent = TRUE)
    if (inherits(time_coef, "try-error") || inherits(ctx_coef, "try-error")) {
      time_coef <- try(tmp2$coefficients[2], silent = TRUE)
      ctx_coef <- try(tmp2$coefficients[3], silent = TRUE)
    }
    tmp2 <- time_coef/2  # Only divide time coefficient by 2 as per theta method
    ctx_effect <- ctx_coef  # Keep context coefficient as is
    context_vectors <- log(ahead::computeattention(x)$context_vectors)
    # Apply modified drift to forecast
    fcast$mean[1] <- fcast$mean[1] + tmp2 * ((1-(1-alpha)^n)/alpha) + ctx_effect*context_vectors
    newx <- c(y, fcast$mean[1])
    for (i in 2:h)
    {
      context_vectors <- log(ahead::computeattention(newx)$context_vectors)
      # Modify drift based on context
      fcast$mean[i] <- fcast$mean[i] + tmp2 * ((i - 1) + (1-(1-alpha)^n)/alpha) + ctx_effect*context_vectors
      newx <- c(newx, fcast$mean[i])
    }
    
    
  } else { # method == "adj"
    context_vectors <- ahead::computeattention(x)$context_vectors
    last_context <- tail(context_vectors, 1)
    # Modify drift based on context
    context_adjusted_drift <- tmp2 * (1 + scale_ctxt * sign(last_context) * 
                                        abs(last_context / mean(abs(y))))
    # Apply modified drift to forecast
    fcast$mean[1] <- fcast$mean[1] + context_adjusted_drift * ((1-(1-alpha)^n)/alpha)
    newx <- c(y, fcast$mean[1])
    for (i in 2:h)
    {
      context_vectors <- ahead::computeattention(newx)$context_vectors
      last_context <- tail(context_vectors, 1)
      # Modify drift based on context
      context_adjusted_drift <- tmp2 * (1 + scale_ctxt * sign(last_context) * 
                                          abs(last_context / mean(abs(newx))))
      fcast$mean[i] <- fcast$mean[i] + context_adjusted_drift * ((i - 1) + (1-(1-alpha)^n)/alpha)
      newx <- c(newx, fcast$mean[i])
    }
  }

  if (seasonal) {
    fcast$mean <- fcast$mean * rep(tail(decomp$seasonal, 
                                        m), trunc(1 + h/m))[1:h]
    fcast$fitted <- fcast$fitted * decomp$seasonal
  }
  fcast$residuals <- origx - fcast$fitted
  fcast.se <- sqrt(fcast$model$sigma2) * sqrt((0:(h - 1)) * 
                                                alpha^2 + 1)
  nconf <- length(level)
  fcast$lower <- fcast$upper <- ts(matrix(NA, nrow = h, ncol = nconf))
  tsp(fcast$lower) <- tsp(fcast$upper) <- tsp(fcast$mean)
  for (i in 1:nconf) {
    zt <- -qnorm(0.5 - level[i]/200)
    fcast$lower[, i] <- fcast$mean - zt * fcast.se
    fcast$upper[, i] <- fcast$mean + zt * fcast.se
  }
  fcast$x <- origx
  fcast$level <- level
  fcast$method <- "Theta"
  fcast$model <- list(alpha = alpha, drift = tmp2, sigma = fcast$model$sigma2)
  fcast$model$call <- match.call()
  return(fcast)
}





