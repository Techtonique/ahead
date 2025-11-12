#' @export
estimate_theta_slope <- function(fit_func,
                                 predict_func,
                                 y,
                                 h = 0,
                                 ...) {
  set.seed(123)
  n <- length(y)
  
  # Extend time range to include forecast horizon
  time_idx <- 0:(n + h - 1)
  n_total <- length(time_idx)
  
  zero <- 1e-4
  eps_factor <- zero^(1 / 3)
  random_covariate <- rnorm(n_total)
  
  # Base data - only first n points have y values
  df <- data.frame(
    y = c(y, rep(NA, h)),
    t = time_idx, 
    z = random_covariate
  )
  
  # Fit model on historical data only (first n rows)
  df_train <- df[1:n, ]
  tmp2 <- try(fit_func(y ~ ., data = df_train, ...), silent = TRUE)
  if (inherits(tmp2, "try-error")) {
    tmp2 <- try(fit_func(x = cbind(df_train$t, df_train$z), 
                         y = df_train$y, ...), silent = TRUE)
    if (inherits(tmp2, "try-error")) {
      tmp2 <- try(fit_func(cbind(df_train$t, df_train$z), 
                           df_train$y, ...), silent = TRUE)
      if (inherits(tmp2, "try-error")) {
        stop("unable to fit the model")
      }
    }
  }
  
  # Adaptive step for ALL time points (historical + future)
  h_eps <- pmax(eps_factor * abs(time_idx), zero)
  double_h <- 2 * h_eps
  
  # Perturbations for all time points
  t_plus  <- time_idx + h_eps
  t_minus <- time_idx - h_eps
  df_plus  <- data.frame(t = t_plus, z = random_covariate)
  df_minus <- data.frame(t = t_minus, z = random_covariate)
  
  # Predict at perturbed points
  fx_plus <- try(as.numeric(predict_func(tmp2, df_plus)), silent = TRUE)
  if (inherits(fx_plus, "try-error"))
    fx_plus <- as.numeric(predict_func(tmp2, as.matrix(df_plus)))
  fx_minus <- try(as.numeric(predict_func(tmp2, df_minus)), silent = TRUE)
  if (inherits(fx_minus, "try-error"))
    fx_minus <- as.numeric(predict_func(tmp2, as.matrix(df_minus)))
  
  # Slopes for ALL time points (historical + future)
  slope <- (fx_plus - fx_minus) / double_h
  
  return(slope)
}

