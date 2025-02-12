#' Conformalized Forecasting using Machine Leaning models
#'
#' @param y A numeric vector or time series of class \code{ts}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags of the input time series considered in the regression
#' @param fit_func Fitting function (Statistical/ML model). Default is Ridge regression.
#' @param predict_func Prediction function (Statistical/ML model)
#' @param coeffs Coefficients of the fitted model. If provided, a linear combination with the coefficients is used to compute the prediction.
#' @param type_pi Type of prediction interval
#' @param B Number of bootstrap replications or number of simulations
#' @param agg "mean" or "median" (aggregation method)
#' @param ... additional parameters passed to the fitting function \code{fit_func}
#'
#' @return
#' @export
#'
#' @examples
#' 
#' res <- ahead::mlf(USAccDeaths, h=10L, lags=15L, type_pi="surrogate", B=250L)
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
mlf <- function(y, h = 5, level = 95, lags = 15L,
                fit_func = ahead::ridge,
                predict_func = predict,
                coeffs = NULL,
                type_pi = c("kde", "surrogate", "blockbootstrap"),
                B = 250L, agg = c("mean", "median"), 
                seed = 123,
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
  y_pred_calibration <- ml_forecast(y = y_train, 
                            h = length(y_calib), 
                            lags=lags, 
                            fit_func = fit_func,
                            predict_func = predict_func,
                            coeffs = coeffs,
                            ...)$mean
  preds_obj <- ml_forecast(y = y_calib, 
                           h = h,  
                           lags = lags, 
                           fit_func = fit_func,
                           predict_func = predict_func,
                           coeffs = coeffs,
                           ...) 
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

    sims <-
      matrix_preds + sd_calibrated_residuals * simulated_scaled_calibrated_residuals

    sims <- ts(sims,
               start = start_preds,
               frequency = frequency(y_train))
    preds_lower <-
      apply(sims, 1, function(x)
        quantile(x, probs = (1 - level / 100) / 2))
    preds_upper <-
      apply(sims, 1, function(x)
        quantile(x, probs = 1 - (1 - level / 100) / 2))

    out <- vector("list", 8) 
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
                        coeffs = NULL, 
                        ...)
{
    df <- as.data.frame(embed(rev(as.numeric(y)), lags + 1L))
    ncol_df <- ncol(df)
    colnames(df) <- c(paste0("lag", rev(seq_len(lags))), "y")    
    if(is.null(coeffs))
    {            
      fit <- try(fit_func(y ~ ., data = df, ...), silent=TRUE)
      if (inherits(fit, "try-error"))
      {
          fit <- fit_func(x = as.matrix(df)[, -ncol_df], 
                          y = df$y, ...)       
      }
      for (i in 1:h)
      {
          newdata <- matrix(df[1, (ncol_df-lags + 1):ncol_df], nrow=1)
          colnames(newdata) <- paste0("lag", rev(seq_len(lags)))
          newdata_df <- data.frame(matrix(as.numeric((newdata)), ncol=lags))
          colnames(newdata_df) <- paste0("lag", rev(seq_len(lags)))
          prediction <- as.numeric(predict(fit, newdata_df)) # /!\ not the same prediction, beware
          df <- data.frame(rbind(as.matrix(cbind(newdata_df, prediction)), as.matrix(df)))
          colnames(df) <- c(paste0("lag", rev(seq_len(lags))), "y")
      }
    } else {
      fit <- list(coefficients = coeffs)
      for (i in 1:h)
      {
          newdata <- matrix(df[1, (ncol_df-lags + 1):ncol_df], nrow=1)
          colnames(newdata) <- paste0("lag", rev(seq_len(lags)))
          newdata_df <- data.frame(matrix(as.numeric((newdata)), ncol=lags))
          colnames(newdata_df) <- paste0("lag", rev(seq_len(lags)))
          prediction <- as.numeric(fit$coefficients %*% as.matrix(newdata_df)) # /!\ not the same prediction, beware
          df <- data.frame(rbind(as.matrix(cbind(newdata_df, prediction)), as.matrix(df)))
          colnames(df) <- c(paste0("lag", rev(seq_len(lags))), "y")
      }
    }    
    return(list(model = fit, mean = head(df$y, h)))
}
