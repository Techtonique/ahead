#' Conformalized Forecasting using Machine Leaning models
#'
#' @param y A numeric vector or time series of class \code{ts}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags of the input time series considered in the regression
#' @param fit_func Fitting function (Statistical/ML model). Default is Ridge regression.
#' @param predict_func Prediction function (Statistical/ML model)
#' @param stack Boolean, use stacking regression or not
#' @param stacking_models A list of \code{fit_func}s and \code{predict_func}s for and ensemble of stacked models (you should set \code{stack=TRUE})
#' @param coeffs Coefficients of the fitted model. If provided, a linear combination with the coefficients is used to compute the prediction.
#' @param type_pi Type of prediction interval
#' @param B Number of bootstrap replications or number of simulations
#' @param agg "mean" or "median" (aggregation method)
#' @param show_progress show progress bar for stacking, if \code{stacking_models} if not \code{NULL}
#' @param ... additional parameters passed to the fitting function \code{fit_func}
#'
#' @return An object of class 'forecast'
#' @export
#'
#' @examples
#' 
#' plot(ahead::mlf(AirPassengers, h=20L))
#'
#' plot(ahead::mlf(AirPassengers, h=25L, lags=20L, fit_func=glmnet::cv.glmnet, stack=TRUE)) 
#'
#' res <- ahead::mlf(USAccDeaths, h=15L, lags=15L, type_pi="surrogate", B=250L)
#' plot(res)
#' 
#' res <- ahead::mlf(USAccDeaths, fit_func = glmnet::cv.glmnet, h=15L, lags=15L, 
#' type_pi="kde", B=250L) 
#' plot(res)
#' 
#' (res <- ahead::mlf(USAccDeaths, fit_func = e1071::svm, h=15L, lags=15L, 
#' type_pi="kde", B=250L)) 
#' plot(res)
#' 
#' res <- ahead::mlf(mdeaths, h=15L, lags=15L, type_pi="surrogate", B=250L)
#' plot(res)
#' 
#' res <- ahead::mlf(fdeaths, fit_func = glmnet::cv.glmnet, h=15L, lags=25L, 
#' type_pi="kde", B=250L) 
#' plot(res)
#' 
#' res <- ahead::mlf(fdeaths, fit_func = randomForest::randomForest, h=15L, lags=15L, 
#' type_pi="kde", B=250L) 
#' plot(res)
#' 
#' 
mlf <- function(y, h = 5, level = 95, lags = 15L,
                fit_func = ahead::ridge,
                predict_func = predict,
                stack = FALSE,
                stacking_models = NULL, 
                coeffs = NULL,
                type_pi = c("surrogate", "bootstrap", "kde"),
                B = 250L, agg = c("mean", "median"), 
                seed = 123,
                show_progress=TRUE,
                ...)
{
  set.seed(seed)
  n <- length(y)
  freq_x <- frequency(y)
  half_n <- ceiling(n/2)
  idx_train <- seq_len(half_n)
  idx_calib <- setdiff(seq_len(n), idx_train)
  splitted_y <- misc::splitts(y)
  y_train <- splitted_y$training
  y_calib <- splitted_y$testing
  type_pi <- match.arg(type_pi)
  agg <- match.arg(agg)
  stacking_results <- list()

  if (!is.null(stacking_models))
  {
    type_pi <- match.arg(type_pi)
    agg <- match.arg(agg)
    n_stacking_models <- length(stacking_models)

    pb <- utils::txtProgressBar(min=0, max=n_stacking_models, style=3)

    j <- 1
    for (i in seq_len(n_stacking_models))
    {
      res <- try(ahead::mlf(y, h=h, lags=lags, 
      fit_func=stacking_models[[i]]$fit_func, 
      predict_func=stacking_models[[i]]$predict_func, 
      stack=stack, stacking_models = NULL, 
      coeffs = NULL, type_pi = type_pi,
      B = B, agg = agg, 
      seed = seed), silent=TRUE)
      if (!inherits(res, "try-error"))
      {
        stacking_results[[j]] <- res 
        j <- j + 1
      } else {
        next 
      }
      utils::setTxtProgressBar(pb, i)
    }   
    close(pb) 

    total_models <- j - 1
    lower_bounds <- sapply(seq_len(total_models), function(i) stacking_results[[i]]$lower)
    upper_bounds <- sapply(seq_len(total_models), function(i) stacking_results[[i]]$upper)
    mean_forecasts <- sapply(seq_len(total_models), function(i) stacking_results[[i]]$mean)
    misc::debug_print(lower_bounds)
    misc::debug_print(upper_bounds)
    misc::debug_print(mean_forecasts)
    if (agg == "median")
    {
      lower_bound <- apply(lower_bounds, 1, median)
      upper_bound <- apply(upper_bounds, 1, median)
    } else {
      lower_bound <- rowMeans(lower_bounds)
      upper_bound <- rowMeans(upper_bounds)
    }


    tspx <- tsp(y)
    start_preds <- tspx[2] + 1 / tspx[3]                         
    freq_x <- frequency(y)
    out <- list() 
    class(out) <- "forecast"
    out$mean <- ts(switch(
      agg,
      median = apply(mean_forecasts, 1, median),
      mean = rowMeans(mean_forecasts)
    ),
    start = start_preds,
    frequency = freq_x)    
    out$lower <- ts(lower_bound,
                    start = start_preds,
                    frequency = freq_x)
    out$upper <- ts(upper_bounds,
                    start = start_preds,
                    frequency = freq_x)
    out$x <- y_calib
    out$level <- level 
    out$method <- "Stacked conformalized ML"
    out$model <- stacking_results
    return(out)
  }

  y_pred_calibration <- ml_forecast(y = y_train, 
                            h = length(y_calib), 
                            lags=lags, 
                            fit_func = fit_func,
                            predict_func = predict_func,
                            coeffs = coeffs,
                            ...)$mean
  if (stack == TRUE)
  {
    preds_obj <- ml_forecast(y = y_calib, 
                           h = h,  
                           lags = lags, 
                           fit_func = fit_func,
                           predict_func = predict_func,
                           xreg = y_pred_calibration, 
                           coeffs = coeffs,
                           ...)  
   } else {
    preds_obj <- ml_forecast(y = y_calib, 
                           h = h,  
                           lags = lags, 
                           fit_func = fit_func,
                           predict_func = predict_func,
                           coeffs = coeffs,
                           ...)  
   }
  
  preds <- preds_obj$mean

  tspx <- tsp(y_calib)
  start_preds <- tspx[2] + 1 / tspx[3]                         
  matrix_preds <- replicate(B, preds)                                                  
  calibrated_raw_residuals <- y_calib - y_pred_calibration
  scaled_calib_resids <- base::scale(calibrated_raw_residuals)
  xm <- attr(scaled_calib_resids, "scaled:center")
  xsd <- attr(scaled_calib_resids, "scaled:scale")
  scaled_calibrated_residuals <- base::scale(calibrated_raw_residuals,
                                             center = TRUE,
                                             scale = TRUE)

  if (type_pi == "kde") {        
        simulated_scaled_calibrated_residuals <-
            rgaussiandens(
            scaled_calibrated_residuals,
            n = h,
            p = B,
            seed = seed
            )
        sd_calibrated_residuals <- sd(calibrated_raw_residuals)
    }

    if (type_pi == "surrogate") {
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
            tseries::surrogate(scaled_calibrated_residuals,
                                    ns =
                                        B)[seq_along(h), ]

        sd_calibrated_residuals <- sd(calibrated_raw_residuals)      
    }

    if (type_pi == "bootstrap") {
      freq_calibrated_raw_residuals <- frequency(calibrated_raw_residuals)
      if (length(calibrated_raw_residuals) <= 2 * freq_calibrated_raw_residuals)
        freq_calibrated_raw_residuals <- 1L
      block_size <-
        ifelse(
          freq_calibrated_raw_residuals > 1,
          2 * freq_calibrated_raw_residuals,
          min(8, floor(
            length(calibrated_raw_residuals) / 2
          ))
        )
      block_size <-
        floor(min(
          max(3L, block_size),
          length(calibrated_raw_residuals) - 1L
        ))
        set.seed(seed)
        simulated_scaled_calibrated_residuals <-
          tseries::tsbootstrap(
            scaled_calibrated_residuals,
            nb =
              B,
            b = floor(block_size),
            type =
              "block"
          )[seq_along(h), ]
        sd_calibrated_residuals <- sd(calibrated_raw_residuals)      
    }

    sims <- matrix_preds 
    sims <- sims + sd_calibrated_residuals * simulated_scaled_calibrated_residuals

    sims <- ts(sims,
               start = start_preds,
               frequency = frequency(y_train))
    preds_lower <-
      apply(sims, 1, function(x)
        quantile(x, probs = (1 - level / 100) / 2))
    preds_upper <-
      apply(sims, 1, function(x)
        quantile(x, probs = 1 - (1 - level / 100) / 2))

    out <- list() 
    class(out) <- "forecast"
    out$mean <- ts(switch(
      agg,
      median = apply(sims, 1, median),
      mean = apply(sims, 1, mean)
    ),
    start = start_preds,
    frequency = freq_x)    
    out$lower <- ts(preds_lower,
                    start = start_preds,
                    frequency = freq_x)
    out$upper <- ts(preds_upper,
                    start = start_preds,
                    frequency = freq_x)
    out$x <- y_calib
    out$level <- level 
    out$method <- "conformalized ML"
    out$model <- preds_obj$model
    out$residuals <- ts(
      calibrated_raw_residuals,
      start = start(y_calib),
      frequency = frequency(y_train)
    ) # /!\ not the same residuals, beware
    out$fitted <- ts(
      y_pred_calibration,
      start = start(y_calib),
      frequency = frequency(y_train)
    ) # /!\ not the same fitted, beware
    out$sims <- sims
    return(out)
}

ml_forecast <- function(y, h, 
                        lags=1, 
                        fit_func = ahead::ridge,
                        predict_func = predict,
                        xreg = NULL, 
                        coeffs = NULL, 
                        ...)
{  
  df <- as.data.frame(embed(rev(as.numeric(y)), lags + 1L))
  ncol_df <- ncol(df)
  colnames(df) <- c(paste0("lag", rev(seq_len(lags))), "y")    
 
  if (!is.null(xreg))
  {
    df_xreg <- as.data.frame(embed(rev(as.numeric(xreg)), lags + 1L))
    ncol_df_xreg <- ncol(df_xreg)
    colnames(df_xreg) <- c(paste0("xreg_lag", rev(seq_len(lags))), "xreg")    
  }  else {
    df_xreg <- NULL 
  }
 
  if (!is.null(xreg))
  {
    fit <- try(fit_func(y ~ ., data = cbind.data.frame(df, df_xreg), ...), silent=TRUE)
    if (inherits(fit, "try-error"))
    {      
      idx_y <- which(colnames(df) == "y")
      fit <- fit_func(x = cbind(as.matrix(df)[, -idx_y], as.matrix(df_xreg)), 
                    y = df$y, ...)             
    }
  } else {
    fit <- try(fit_func(y ~ ., data = df, ...), silent=TRUE)
    if (inherits(fit, "try-error"))
    {      
      idx_y <- which(colnames(df) == "y")
      fit <- fit_func(x = as.matrix(df)[, -idx_y], 
                      y = df$y, ...)             
    }
  }
  
  forecasts <- numeric(h)
  
  # Convert to matrix once for faster operations
  df_matrix <- as.matrix(df)
  col_names <- colnames(df)

  if (!is.null(xreg))
  {
    df_xreg_matrix <- as.matrix(df_xreg)
    col_names_xreg <- colnames(df_xreg)
  }

  if (is.null(xreg))
  {
    for (i in 1:h)
    {
      # Extract newdata exactly as original
      newdata <- matrix(df_matrix[1, (ncol_df-lags + 1):ncol_df], nrow=1)
      colnames(newdata) <- paste0("lag", rev(seq_len(lags)))
      newdata_df <- data.frame(matrix(as.numeric((newdata)), ncol=lags))
      colnames(newdata_df) <- paste0("lag", rev(seq_len(lags)))
      prediction <- as.numeric(predict_func(fit, newdata_df))
      forecasts[i] <- prediction
      # Create new row and rbind to matrix (faster than data.frame rbind)
      new_row <- c(as.numeric(newdata_df), prediction)
      df_matrix <- rbind(new_row, df_matrix)
    } 
  } else {
    for (i in 1:h)
    {
      # Extract newdata exactly as original      
      newdata_y <- matrix(df_matrix[1, (ncol_df-lags + 1):ncol_df], nrow=1)
      newdata <- cbind(newdata_y, 
                       matrix(df_xreg_matrix[1, (ncol_df-lags):ncol_df], nrow=1))
      colnames(newdata) <- c(paste0("lag", rev(seq_len(lags))), c(paste0("xreg_lag", rev(seq_len(lags))), "xreg")) 
      newdata_y_df <- data.frame(matrix(as.numeric((newdata_y)), ncol=lags))    
      colnames(newdata_y_df) <- paste0("lag", rev(seq_len(lags)))
      newdata_df <- data.frame(matrix(as.numeric((newdata)), ncol=2*lags + 1))
      colnames(newdata_df) <- c(paste0("lag", rev(seq_len(lags))), c(paste0("xreg_lag", rev(seq_len(lags))), "xreg"))
      prediction <- as.numeric(predict_func(fit, newdata_df))
      forecasts[i] <- prediction
      # Create new row and rbind to matrix (faster than data.frame rbind)
      new_row <- c(as.numeric(newdata_y_df), prediction)
      df_matrix <- rbind(new_row, df_matrix)
    }
  }   
  
  return(list(model = fit, mean = forecasts))
}