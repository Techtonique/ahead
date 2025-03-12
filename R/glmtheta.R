#' @title Generalized Linear Model Theta Forecast
#' @description This function implements the Theta method using a Generalized Linear Model (GLM)
#' @param y The time series data
#' @param h The number of periods to forecast
#' @param level The confidence level for the forecast intervals
#' @param fit_func The function to use for fitting the GLM
#' @param fan Logical flag for fan plot
#' @param x The time series data
#' @param attention Logical flag for using attention mechanism
#' @param ... Additional arguments to pass to the fit_func
#' @return A forecast object
#' @export
glmthetaf <- function (y, h = ifelse(frequency(y) > 1, 2 * frequency(y), 10), 
                       level = c(80, 95), fit_func = lsfit, fan = FALSE, x = y, 
                       attention = FALSE, method=c("base", "adj"),
                       ...) 
{
  method <- match.arg(method) 
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
  } else { # attention is TRUE
      if (method == "base")
      {
        df <- cbind.data.frame(df, ctx = ahead::computeattention(x)$context_vectors)
        tmp2 <- try(fit_func(y ~ ., data = df, ...), silent = FALSE)
        if (inherits(tmp2, "try-error")) {
          X <- cbind.data.frame(time_idx, df$ctx)
          tmp2 <- try(fit_func(x = as.matrix(X), y = as.numeric(x), ...), silent = FALSE)
          if (inherits(tmp2, "try-error")) {
            stop("Unable to fit linear trend with attention")
          }
        }
      } else { # method == "adj": same as attention == FALSE to avoid adjusting the trend twice
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
    time_coef <- try(coef(tmp2)[2], silent = TRUE)
    ctx_coef <- try(coef(tmp2)[3], silent = TRUE)
    if (inherits(time_coef, "try-error") || inherits(ctx_coef, "try-error")) {
      time_coef <- try(tmp2$coefficients[2], silent = TRUE)
      ctx_coef <- try(tmp2$coefficients[3], silent = TRUE)
    }
    tmp2 <- time_coef/2  # Only divide time coefficient by 2 as per theta method
    ctx_effect <- ctx_coef  # Keep context coefficient as is
    fcast$mean <- fcast$mean + tmp2 * (0:(h-1) + (1-(1-alpha)^n)/alpha)
  } else { # method == "adj"
    context_vectors <- ahead::computeattention(x)$context_vectors
    last_context <- tail(context_vectors, 1)
    # Modify drift based on context
    context_adjusted_drift <- tmp2 * (1 + sign(last_context) * 
                                        abs(last_context / mean(abs(y))))
    # Apply modified drift to forecast
    fcast$mean[1] <- fcast$mean[1] + context_adjusted_drift * ((1-(1-alpha)^n)/alpha)
    newx <- c(y, fcast$mean[1])
    for (i in 2:h)
    {
      context_vectors <- ahead::computeattention(newx)$context_vectors
      last_context <- tail(context_vectors, 1)
      # Modify drift based on context
      context_adjusted_drift <- tmp2 * (1 + sign(last_context) * 
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





