  # Simple 2-Level Top-Down Hierarchical Forecasting with Sequential Split Conformal
  # Returns simulations at both total and bottom levels

  #' Create simple 2-level hierarchical time series
  #' @param bts Bottom-level time series matrix (T x m)
  #' @return List with total and bottom series
  #' @export 
  create_2level_hts <- function(bts) {
    if (!is.matrix(bts)) bts <- as.matrix(bts)
    
    list(
      bts = bts,
      total = rowSums(bts),
      proportions = colMeans(bts / rowSums(bts), na.rm = TRUE)
    )
  }

  #' Simple ETS forecast function
  #' @param y Time series vector
  #' @param h Forecast horizon
  #' @return Forecast vector
  simple_forecast <- function(y, h) {
    tryCatch({
      fit <- ets(ts(y))
      as.numeric(forecast(fit, h = h)$mean)
    }, error = function(e) {
      rep(tail(y, 1), h)  # Fallback to naive
    })
  }

  #' Top-down forecast using historical proportions
  #' @param hts 2-level HTS object
  #' @param h Forecast horizon
  #' @return List with total and bottom forecasts
  #' @export 
  topdown_forecast <- function(hts, h) {
    # Forecast total
    total_fc <- simple_forecast(hts$total, h)
    
    # Disaggregate using historical proportions
    bottom_fc <- matrix(total_fc, nrow = h, ncol = ncol(hts$bts)) * 
      matrix(hts$proportions, nrow = h, ncol = ncol(hts$bts), byrow = TRUE)
    
    list(total = total_fc, bottom = bottom_fc)
  }

  #' Sequential split conformal prediction for hierarchical forecasting
  #' Returns simulations at both total and bottom levels
  #' @param hts 2-level HTS object
  #' @param h Forecast horizon
  #' @param split_ratio Training/calibration split ratio
  #' @param n_sim Number of simulation samples
  #' @param alpha Significance level for prediction intervals
  #' @return Forecast with prediction intervals and simulations
  #' @export 
  sequential_conformal_hts <- function(hts, h = 12, split_ratio = 0.6, 
                                      n_sim = 1000, alpha = 0.1) {  
    n_obs <- nrow(hts$bts)
    n_series <- ncol(hts$bts)  
    # STEP 1: Sequential split
    split_point <- floor(n_obs * split_ratio)
    split_point <- max(h, min(split_point, n_obs - h))  
    # Training data (first part)
    train_bts <- hts$bts[1:split_point, , drop = FALSE]
    train_hts <- create_2level_hts(train_bts)  
    # Calibration data (second part)
    cal_start <- split_point + 1
    cal_end <- min(split_point + h, n_obs)
    cal_length <- cal_end - cal_start + 1  
    cal_bts_true <- hts$bts[cal_start:cal_end, , drop = FALSE]
    cal_total_true <- rowSums(cal_bts_true)  
    # STEP 2: Forecast calibration period using training data
    cal_forecast <- topdown_forecast(train_hts, cal_length)  
    # STEP 3: Calculate calibration residuals
    cal_residuals_bottom <- cal_bts_true - cal_forecast$bottom[1:cal_length, ]
    cal_residuals_total <- cal_total_true - cal_forecast$total[1:cal_length]  
    # STEP 4: Train final model on calibration data
    final_bts <- hts$bts[cal_start:n_obs, , drop = FALSE]
    final_hts <- create_2level_hts(final_bts)
    final_forecast <- topdown_forecast(final_hts, h)  
    # STEP 5: Generate simulations
    sim_bottom <- array(NA, dim = c(h, n_series, n_sim))
    sim_total <- matrix(NA, nrow = h, ncol = n_sim)  
    for (i in 1:n_sim) {
      # Sample residuals with replacement for each series
      boot_residuals <- matrix(NA, nrow = h, ncol = n_series)
      
      for (j in 1:n_series) {
        series_residuals <- cal_residuals_bottom[, j]
        series_residuals <- series_residuals[!is.na(series_residuals)]
        
        if (length(series_residuals) > 0) {
          # Sample residuals with replacement
          boot_residuals[, j] <- sample(series_residuals, h, replace = TRUE)
        } else {
          boot_residuals[, j] <- 0
        }
      }    
      # Add to forecasts (with random sign for symmetry)
      signs <- matrix(sample(c(-1, 1), h * n_series, replace = TRUE), 
                      nrow = h, ncol = n_series)    
      # Create simulation
      sim_bottom[, , i] <- final_forecast$bottom + signs * boot_residuals    
      # Ensure non-negative values (if appropriate for your data)
      sim_bottom[, , i] <- pmax(sim_bottom[, , i], 0)
      
      # Calculate total from bottom simulations (maintains hierarchy)
      sim_total[, i] <- rowSums(sim_bottom[, , i])
    }
    
    # STEP 6: Calculate prediction intervals
    lower_q <- alpha / 2
    upper_q <- 1 - alpha / 2
    
    # For bottom series
    lower_bottom <- apply(sim_bottom, c(1, 2), quantile, probs = lower_q, na.rm = TRUE)
    upper_bottom <- apply(sim_bottom, c(1, 2), quantile, probs = upper_q, na.rm = TRUE)
    
    # For total series
    lower_total <- apply(sim_total, 1, quantile, probs = lower_q, na.rm = TRUE)
    upper_total <- apply(sim_total, 1, quantile, probs = upper_q, na.rm = TRUE)
    
    # Calculate mean forecasts from simulations
    mean_bottom <- apply(sim_bottom, c(1, 2), mean, na.rm = TRUE)
    mean_total <- apply(sim_total, 1, mean, na.rm = TRUE)
    
    # Return comprehensive results with simulations
    list(
      point_forecasts = list(
        total = final_forecast$total,
        bottom = final_forecast$bottom
      ),
      mean_forecasts = list(
        total = mean_total,
        bottom = mean_bottom
      ),
      prediction_intervals = list(
        total = cbind(lower = lower_total, upper = upper_total),
        bottom = list(lower = lower_bottom, upper = upper_bottom)
      ),
      simulations = list(
        total = sim_total,
        bottom = sim_bottom
      ),
      calibration_residuals = list(
        bottom = cal_residuals_bottom,
        total = cal_residuals_total
      ),
      split_info = list(
        split_point = split_point,
        calibration_length = cal_length,
        split_ratio = split_ratio,
        n_sim = n_sim,
        alpha = alpha
      )
    )
  }

  #' Plot forecast results with simulation intervals
  plot_hts_forecast <- function(result, hts, series_type = "total", series_idx = 1) {
    
    if (series_type == "total") {
      historical <- hts$total
      point_fc <- result$point_forecasts$total
      mean_fc <- result$mean_forecasts$total
      lower <- result$prediction_intervals$total[, "lower"]
      upper <- result$prediction_intervals$total[, "upper"]
      title <- "Total Series Forecast"
    } else {
      historical <- hts$bts[, series_idx]
      point_fc <- result$point_forecasts$bottom[, series_idx]
      mean_fc <- result$mean_forecasts$bottom[, series_idx]
      lower <- result$prediction_intervals$bottom$lower[, series_idx]
      upper <- result$prediction_intervals$bottom$upper[, series_idx]
      title <- paste("Bottom Series", series_idx, "Forecast")
    }
    
    h <- length(point_fc)
    n_hist <- length(historical)
    
    plot(1:n_hist, historical, type = "l", lwd = 2, col = "black",
        xlim = c(n_hist - 50, n_hist + h),
        ylim = range(c(tail(historical, 50), point_fc, lower, upper), na.rm = TRUE),
        xlab = "Time", ylab = "Value", main = title)
    
    # Add point forecast
    lines((n_hist + 1):(n_hist + h), point_fc, col = "blue", lwd = 2)
    
    # Add mean forecast from simulations
    lines((n_hist + 1):(n_hist + h), mean_fc, col = "red", lwd = 2, lty = 2)
    
    # Add prediction intervals
    polygon(c((n_hist + 1):(n_hist + h), rev((n_hist + 1):(n_hist + h))),
            c(lower, rev(upper)), col = rgb(0.5, 0.5, 1, 0.3), border = NA)
    
    # Add split line
    abline(v = result$split_info$split_point + 0.5, col = "red", lty = 2)
    
    legend("topleft", 
          legend = c("Historical", "Point Forecast", "Mean Forecast", 
                      paste0((1-result$split_info$alpha)*100, "% PI"), "Train/Test Split"),
          col = c("black", "blue", "red", "lightblue", "red"),
          lty = c(1, 1, 2, 1, 2), lwd = c(2, 2, 2, 8, 1))
  }

  #' Plot simulation paths
  #' @export 
  plot_simulations <- function(result, series_type = "total", series_idx = 1, 
                              n_paths = 50, main = NULL) {
    
    if (series_type == "total") {
      sim_data <- result$simulations$total
      title <- if (is.null(main)) "Total Series Simulation Paths" else main
      ylab <- "Total Value"
    } else {
      sim_data <- result$simulations$bottom[, series_idx, ]
      title <- if (is.null(main)) paste("Bottom Series", series_idx, "Simulation Paths") else main
      ylab <- paste("Series", series_idx, "Value")
    }
    
    h <- nrow(sim_data)
    n_sim <- ncol(sim_data)
    
    # Plot first simulation path
    plot(1:h, sim_data[, 1], type = "l", col = rgb(0.5, 0.5, 0.5, 0.3),
        ylim = range(sim_data, na.rm = TRUE),
        xlab = "Forecast Horizon", ylab = ylab, main = title)
    
    # Add additional paths
    n_plot <- min(n_paths, n_sim)
    for (i in 2:n_plot) {
      lines(1:h, sim_data[, i], col = rgb(0.5, 0.5, 0.5, 0.3))
    }
    
    # Add mean path
    lines(1:h, rowMeans(sim_data), col = "red", lwd = 2)
    
    # Add prediction intervals
    lower <- apply(sim_data, 1, quantile, probs = result$split_info$alpha/2, na.rm = TRUE)
    upper <- apply(sim_data, 1, quantile, probs = 1 - result$split_info$alpha/2, na.rm = TRUE)
    
    polygon(c(1:h, rev(1:h)), c(lower, rev(upper)), 
            col = rgb(1, 0, 0, 0.2), border = NA)
    
    legend("topleft", 
          legend = c("Simulation Paths", "Mean Path", 
                      paste0((1-result$split_info$alpha)*100, "% PI")),
          col = c("gray", "red", "pink"),
          lty = c(1, 1, 1), lwd = c(1, 2, 8))
  }

  #' Example usage
  example_2level_forecast <- function() {
    # Create example data
    set.seed(123)
    n_obs <- 120
    n_series <- 3
    
    # Generate correlated bottom series
    bts <- matrix(NA, nrow = n_obs, ncol = n_series)
    base <- 50 + cumsum(rnorm(n_obs, 0, 1))
    for (i in 1:n_series) {
      bts[, i] <- base * (0.8 + 0.1 * i) + rnorm(n_obs, 0, 2)
    }
    
    # Create HTS object
    hts <- create_2level_hts(bts)
    
    # Run sequential conformal forecast
    result <- sequential_conformal_hts(hts, h = 12, n_sim = 1000)
    
    # Plot results
    par(mfrow = c(2, 2))
    plot_hts_forecast(result, hts, "total")
    for (i in 1:min(2, n_series)) {
      plot_hts_forecast(result, hts, "bottom", i)
    }
    par(mfrow = c(1, 1))
    
    # Plot simulation paths
    par(mfrow = c(2, 2))
    plot_simulations(result, "total")
    for (i in 1:min(2, n_series)) {
      plot_simulations(result, "bottom", i)
    }
    par(mfrow = c(1, 1))
    
    # Display simulation statistics
    cat("=== Simulation Statistics ===\n")
    cat("Total simulations:", result$split_info$n_sim, "\n")
    cat("Total forecast mean:", round(mean(result$simulations$total), 2), "\n")
    cat("Total forecast SD:", round(sd(result$simulations$total), 2), "\n")
    
    return(result)
  }

  # Function to extract simulation data for further analysis
  get_simulation_data <- function(result, series_type = "total", series_idx = 1) {
    if (series_type == "total") {
      return(result$simulations$total)
    } else {
      return(result$simulations$bottom[, series_idx, ])
    }
  }
