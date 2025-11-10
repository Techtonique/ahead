estimate_theta_slope <- function(fit_func,
                                 predict_func,
                                 y,
                                 ...) {
  set.seed(123)
  n <- length(y)
  time_idx <- 0:(n - 1)
  zero <- 1e-4
  eps_factor <- zero^(1 / 3)
  random_covariate <- rnorm(n)
  # Base data
  df <- data.frame(y = y, t = time_idx, 
                   z = random_covariate)
  # Fit model (try formula interface, then x/y)
  tmp2 <- try(fit_func(y ~ ., data = df, ...), silent = FALSE)
  if (inherits(tmp2, "try-error")) {
    tmp2 <- try(fit_func(x = cbind(time_idx, random_covariate), 
                         y = y, ...), silent = FALSE)
    if (inherits(tmp2, "try-error"))
    {
      tmp2 <- try(fit_func(cbind(time_idx, random_covariate), 
                           y, ...), silent = FALSE)
      if (inherits(tmp2, "try-error")) {
        stop("unable to fit the model, check the 'estimate_theta_slope' function's code")
      }
    }
  }
  # Adaptive step
  h_eps <- pmax(eps_factor * abs(time_idx), zero)
  double_h <- 2 * h_eps
  # Perturbations
  t_plus  <- time_idx + h_eps
  t_minus <- time_idx - h_eps
  df_plus  <- data.frame(t = t_plus, z = random_covariate)
  df_minus <- data.frame(t = t_minus, z = random_covariate)
  # Robust predictions
  fx_plus <- try(as.numeric(predict_func(tmp2, df_plus)), silent = FALSE)
  if (inherits(fx_plus, "try-error"))
    fx_plus <- as.numeric(predict_func(tmp2, as.matrix(df_plus)))
  fx_minus <- try(as.numeric(predict_func(tmp2, df_minus)), silent = FALSE)
  if (inherits(fx_minus, "try-error"))
    fx_minus <- as.numeric(predict_func(tmp2, as.matrix(df_minus)))
  slope <- (fx_plus - fx_minus) / double_h
  return(median(slope))
}
