#' Ridge2 model
#'
#' Random Vector functional link network model with 2 regularization parameters
#'
#' @param y A univariate of multivariate time series of class \code{ts} (preferred) or a \code{matrix}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags
#' @param xreg External regressors. A data.frame (preferred) or a \code{matrix}
#' @param nb_hidden Number of nodes in hidden layer
#' @param nodes_sim Type of simulation for nodes in the hidden layer
#' @param activ Activation function
#' @param a Hyperparameter for activation function "leakyrelu", "elu"
#' @param lambda_1 Regularization parameter for original predictors
#' @param lambda_2 Regularization parameter for transformed predictors
#' @param dropout dropout regularization parameter (dropping nodes in hidden layer)
#' @param seed Reproducibility seed for `nodes_sim == unif`
#' @param type_forecast Recursive or direct forecast
#' @param type_pi Type of prediction interval currently "gaussian", "bootstrap",
#' "blockbootstrap", "movingblockbootstrap", "rvinecopula" (with Gaussian and empirical margins for now)
#' @param block_length Length of block for circular or moving block bootstrap
#' @param margins Distribution of margins: "gaussian", "empirical" for \code{type_pi == "rvinecopula"}
#' @param seed Reproducibility seed for random stuff
#' @param B Number of bootstrap replications or number of simulations (yes, 'B' is unfortunate)
#' @param type_aggregation Type of aggregation, ONLY for bootstrapping; either "mean" or "median"
#' @param centers Number of clusters for \code{type_clustering}
#' @param type_clustering "kmeans" (K-Means clustering) or "hclust" (Hierarchical clustering)
#' @param ym Univariate time series (\code{stats::ts}) of yield to maturities with
#' \code{frequency = frequency(y)} and \code{start = tsp(y)[2] + 1 / frequency(y)}.
#' Default is \code{NULL}.
#' @param cl An integer; the number of clusters for parallel execution, for bootstrap
#' @param show_progress A boolean; show progress bar for bootstrapping? Default is TRUE.
#' @param ... Additional parameters to be passed to \code{\link{kmeans}} or \code{\link{hclust}}
#'
#' @return An object of class "mtsforecast"; a list containing the following elements:
#'
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{sims}{Model simulations for bootstrapping (basic, or block)}
#' \item{x}{The original time series}
#' \item{residuals}{Residuals from the fitted model}
#' \item{coefficients}{Regression coefficients for \code{type_pi == 'gaussian'} for now}
#'
#' @author T. Moudiki
#'
#' @references
#'
#' Moudiki, T., Planchet, F., & Cousin, A. (2018).
#' Multiple time series forecasting using quasi-randomized
#' functional link neural networks. Risks, 6(1), 22. \cr
#'
#' @export
#'
#' @examples
#'
#' require(fpp)
#'
#' print(ahead::ridge2f(fpp::insurance)$mean)
#' print(ahead::ridge2f(fpp::usconsumption)$lower)
#'
#' res <- ahead::ridge2f(fpp::insurance, h=10, lags=2)
#' par(mfrow=c(1, 2))
#' plot(res, "Quotes")
#' plot(res, "TV.advert")
#'
#' # include a trend (just for the example)
#' xreg <- as.numeric(time(fpp::insurance))
#' res2 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
#' h=10, lags=2)
#' par(mfrow=c(1, 2))
#' plot(res2, "Quotes")
#' plot(res2, "TV.advert")
#'
#' # block bootstrap
#' xreg <- as.numeric(time(fpp::insurance))
#' res3 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
#'                       h=10, lags=1L, type_pi = "bootstrap", B=10)
#' res5 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
#'                       h=10, lags=1L, type_pi = "blockbootstrap", B=10,
#'                       block_length = 4)
#'
#' print(res3$sims[[2]])
#' print(res5$sims[[2]])
#'
#' par(mfrow=c(2, 2))
#' plot(res3, "Quotes")
#' plot(res3, "TV.advert")
#' plot(res5, "Quotes")
#' plot(res5, "TV.advert")
#'
#'
#' res4 <- ahead::ridge2f(fpp::usconsumption, h=20, lags=2L,
#' lambda_2=1)
#' par(mfrow=c(1, 2))
#' plot(res4, "income")
#' plot(res4, "consumption")
#'
#'
#' # moving block bootstrap
#' xreg <- as.numeric(time(fpp::insurance))
#' res6 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
#'                       h=10, lags=1L,
#'                       type_pi = "movingblockbootstrap", B=10,
#'                       block_length = 4)
#'
#' print(res6$sims[[2]])
#'
#' par(mfrow=c(1, 2))
#' plot(res6, "Quotes")
#' plot(res6, "TV.advert")
#'
#'
ridge2f <- function(y,
                    h = 5,
                    level = 95,
                    xreg = NULL,
                    lags = 1,
                    nb_hidden = 5,
                    nodes_sim = c("sobol", "halton", "unif"),
                    activ = c("relu", "sigmoid", "tanh",
                              "leakyrelu", "elu", "linear"),
                    a = 0.01,
                    lambda_1 = 0.1,
                    lambda_2 = 0.1,
                    dropout = 0,
                    type_forecast = c("recursive", "direct"),
                    type_pi = c(
                      "gaussian",
                      "bootstrap",
                      "blockbootstrap",
                      "movingblockbootstrap",
                      "rvinecopula"
                    ),
                    block_length = NULL,
                    margins = c("gaussian", "empirical"),
                    seed = 1,
                    B = 100L,
                    type_aggregation = c("mean", "median"),
                    centers = NULL,
                    type_clustering = c("kmeans", "hclust"),
                    ym = NULL,
                    cl = 1L,
                    show_progress = TRUE,
                    ...)
{
  if(is_package_available("randtoolbox") == FALSE)
    install.packages("randtoolbox",
                     repos = c(CRAN = "https://cloud.r-project.org"))

  if(is_package_available("VineCopula") == FALSE)
    install.packages("VineCopula",
                     repos = c(CRAN = "https://cloud.r-project.org"))

  if (is.null(dim(y)))
  {
    dimensionality <- "univariate"

    is_ts_y_uni <- is.ts(y)

    if (is_ts_y_uni)
    {
      start_y_uni <- start(y)
      freq_y_uni <- frequency(y)
    }
    y <- cbind(y, seq_along(y))
    colnames(y) <- c("y", "trend_univariate")
    if (is_ts_y_uni){
      y <- ts(y, start = start_y_uni,
              frequency = freq_y_uni)
    }
  } else {
    dimensionality <- "multivariate"
  }

  stopifnot(floor(lags) == lags)

  n_series <- ncol(y)

  colnames_y <- colnames(y)
  if (is.null(colnames_y))
  {
    colnames_y <- colnames(y) <- paste0("y", 1:n_series)
  }

  if (!is.ts(y))
  {
    y <- ts(y)
  }

  use_xreg <- FALSE
  use_clustering <- FALSE

  if (!is.null(xreg))
  {
    if (is.null(ncol(xreg)))
      xreg <- as.matrix(xreg)

    stopifnot(identical(nrow(xreg), nrow(y)))

    use_xreg <- TRUE

    n_xreg <- dim(xreg)[2]
    colnames_xreg <- colnames(xreg)
    if (is.null(colnames_xreg))
    {
      colnames(xreg) <- paste0("xreg_", 1:n_xreg)
      colnames_xreg <- colnames(xreg)
    } else {
      colnames(xreg) <- paste0("xreg_", colnames_xreg)
      colnames_xreg <- colnames(xreg)
    }
    y <- cbind(y, xreg)
    colnames(y) <- c(colnames_y, colnames_xreg)
  }

  if (!is.null(centers))
  {
    use_clustering <- TRUE
    centers <- floor(max(2L, min(centers, nrow(y) - 1L)))
    matrix_clusters <- get_clusters(y,
                                    centers = centers,
                                    type_clustering = match.arg(type_clustering),
                                    seed = seed,
                                    ...) # return a matrix, not a ts # additional params for clustering
    colnames_clustering <- colnames(matrix_clusters)
    colnames_y_clustering <- colnames(y)
    y <- cbind(y, matrix_clusters)
    colnames(y) <- c(colnames_y_clustering, colnames_clustering)
  }

  series_names <- colnames(y)
  type_pi <- match.arg(type_pi)
  freq_x <- frequency(y)
  start_x <- start(y)
  start_preds <- tsp(y)[2] + 1 / freq_x

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_forecast <- match.arg(type_forecast)

  # Fitting a regularized regression  model to multiple time series
  if (!identical(type_pi, "splitconformal"))
  {
    fit_obj <- fit_ridge2_mts(
      y,
      lags = lags,
      nb_hidden = nb_hidden,
      nodes_sim = nodes_sim,
      activ = activ,
      a = a,
      lambda_1 = lambda_1,
      lambda_2 = lambda_2,
      dropout = dropout,
      centers = centers,
      type_clustering = type_clustering,
      seed = seed,
      type_pi = type_pi,
      margins = margins
    )
  } else { # if (type_pi == "splitconformal") # experimental
    # Split the training data
    y_train_calibration <- splitts(y, split_prob=0.5)
    y_train <- y_train_calibration$training 
    y_calibration <- y_train_calibration$testing
    h_calibration <- nrow(y_calibration)

    # First fit on training data
    fit_obj_train <- fit_ridge2_mts(
      y_train,
      lags = lags,
      nb_hidden = nb_hidden,
      nodes_sim = nodes_sim,
      activ = activ,
      a = a,
      lambda_1 = lambda_1,
      lambda_2 = lambda_2,
      dropout = dropout,
      centers = centers,
      type_clustering = type_clustering,
      seed = seed
    )

    # Get predictions on calibration set
    y_pred_calibration <- ridge2f(
      y_train,
      h = h_calibration,
      lags = lags,
      nb_hidden = nb_hidden,
      nodes_sim = nodes_sim,
      activ = activ,
      a = a,
      lambda_1 = lambda_1,
      lambda_2 = lambda_2,
      dropout = dropout,
      centers = centers,
      type_clustering = type_clustering,
      seed = seed
    )$mean

    # Calculate residuals on calibration set
    matrix_y_calibration <- matrix(as.numeric(y_calibration), ncol = ncol(y))
    matrix_y_pred_calibration <- matrix(as.numeric(y_pred_calibration), ncol = ncol(y))
    abs_residuals <- abs(matrix_y_calibration - matrix_y_pred_calibration)
    
    # Get conformal quantile for each series
    quantile_absolute_residuals_conformal <- apply(abs_residuals, 2, function(x) {
      quantile_scp(
        abs_residuals = x,
        alpha = (1 - level/100)
      )
    })

    # Final fit and forecast on full calibration set
    preds <- ridge2f(
      y_calibration, 
      h = h,
      lags = lags,
      nb_hidden = nb_hidden,
      nodes_sim = nodes_sim,
      activ = activ,
      a = a,
      lambda_1 = lambda_1,
      lambda_2 = lambda_2,
      dropout = dropout,
      centers = centers,
      type_clustering = type_clustering,
      seed = seed
    )$mean

    # Create output with conformal intervals (now using matrix operations)
    preds_matrix <- matrix(as.numeric(preds), ncol = ncol(y))
    lower_matrix <- sweep(preds_matrix, 2, quantile_absolute_residuals_conformal, "-")
    upper_matrix <- sweep(preds_matrix, 2, quantile_absolute_residuals_conformal, "+")

    out <- list(
      mean = preds,
      lower = ts(lower_matrix, start = tsp(preds)[1], frequency = tsp(preds)[3]),
      upper = ts(upper_matrix, start = tsp(preds)[1], frequency = tsp(preds)[3]),
      sims = NULL,
      x = y,
      level = level,
      method = "ridge2",
      residuals = ts(abs_residuals, start = start_x),
      coefficients = fit_obj_train$coef,
      loocv = fit_obj_train$loocv,
      weighted_loocv = fit_obj_train$weighted_loocv, 
      loocv_per_series = fit_obj_train$loocv_per_series
    )

    # Handle univariate case
    if (dimensionality == "univariate") {
      for (i in 1:length(out)) {
        try_delete_trend <- try(delete_columns(out[[i]], "trend_univariate"), silent = TRUE)
        if (!inherits(try_delete_trend, "try-error") && !is.null(out[[i]])) {
          out[[i]] <- try_delete_trend
        }
      }
      return(structure(out, class = "forecast"))
    }

    return(structure(out, class = "mtsforecast"))
  }

  if (identical(type_pi, "gaussian"))
  {
    preds <- ts(
      data = fcast_ridge2_mts(
        fit_obj,
        h = h,
        type_forecast = type_forecast,
        level = level
      ),
      start = start_preds,
      frequency = freq_x
    )

    # strong Gaussian assumption
    sds <- apply(fit_obj$resids, 2, sd)
    qtsd <- -qnorm(0.5 - level / 200) * sds
    qtsd <- t(replicate(nrow(preds), qtsd))

    # Forecast from fit_obj
    out <- list(
      mean = preds,
      lower = preds - qtsd,
      upper = preds + qtsd,
      sims = NULL,
      x = y,
      level = level,
      method = "ridge2",
      residuals = ts(fit_obj$resids,
                     start = start_x),
      coefficients = fit_obj$coef,
      loocv = fit_obj$loocv,
      weighted_loocv = fit_obj$weighted_loocv,
      loocv_per_series = fit_obj$loocv_per_series
    )

    if (use_xreg || use_clustering)
    {
      for (i in 1:length(out))
      {
        try_delete_xreg <-
          try(delete_columns(out[[i]], "xreg_"), silent = TRUE)
        if (!inherits(try_delete_xreg, "try-error") &&
            !is.null(out[[i]]))
        {
          out[[i]] <- try_delete_xreg
        }
      }
    }

    if (dimensionality == "univariate")
    {
      for (i in 1:length(out))
      {
        try_delete_trend <-
          try(delete_columns(out[[i]], "trend_univariate"), silent = TRUE)
        if (!inherits(try_delete_trend, "try-error") &&
            !is.null(out[[i]]))
        {
          out[[i]] <- try_delete_trend
        }
      }
      return(structure(out, class = "forecast"))
    }

    return(structure(out, class = "mtsforecast"))
  }

  if (type_pi %in% c("bootstrap", "blockbootstrap",
                     "movingblockbootstrap"))
  {
    cl <- floor(min(max(cl, 0L), parallel::detectCores()))

    if (cl <= 1L)
    {

      if (show_progress == FALSE)
      {

        sims <- lapply(1:B,
                       function(i)
                         ts(
                           fcast_ridge2_mts(
                             fit_obj,
                             h = h,
                             type_forecast = type_forecast,
                             bootstrap = TRUE,
                             type_bootstrap = type_pi,
                             block_length = block_length,
                             seed = seed + i * 100
                           ),
                           start = start_preds,
                           frequency = freq_x
                         ))

      } else { # show_progress == TRUE

        sims <- base::vector("list", B)

        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)

        for (i in 1:B)
        {
          sims[[i]] <- ts(
            fcast_ridge2_mts(
              fit_obj,
              h = h,
              type_forecast = type_forecast,
              bootstrap = TRUE,
              type_bootstrap = type_pi,
              block_length = block_length,
              seed = seed + i * 100
            ),
            start = start_preds,
            frequency = freq_x
          )
          utils::setTxtProgressBar(pb, i)
        }

        base::close(pb)

      }

    } else { # parallel exec.

      if (show_progress == FALSE)
      {
        cluster <- parallel::makeCluster(getOption("cl.cores", cl))
        sims <-
          parallel::parLapply(
            cl = cluster,
            X = 1:B,
            fun = function(i)
              ts(
                fcast_ridge2_mts(
                  fit_obj,
                  h = h,
                  type_forecast = type_forecast,
                  bootstrap = TRUE,
                  type_bootstrap = type_pi,
                  block_length = block_length,
                  seed = seed + i * 100
                ),
                start = start_preds,
                frequency = freq_x
              )
          )
        parallel::stopCluster(cluster)
      } else { # show_progress == TRUE

        cl_SOCK <- parallel::makeCluster(cl, type = "SOCK")

        doSNOW::registerDoSNOW(cl_SOCK)

        `%op%` <-  foreach::`%dopar%`

        pb <- utils::txtProgressBar(min = 0,
                             max = B,
                             style = 3)

        progress <- function(i) utils::setTxtProgressBar(pb, i)

        opts <- list(progress = progress)

        i <- NULL

        sims <- foreach::foreach(
          i = 1:B,
          #.packages = packages,
          #.combine = rbind,
          .errorhandling = "stop",
          .options.snow = opts,
          .verbose = FALSE
        ) %op% {

          ts(
            fcast_ridge2_mts(
              fit_obj,
              h = h,
              type_forecast = type_forecast,
              bootstrap = TRUE,
              type_bootstrap = type_pi,
              block_length = block_length,
              seed = seed + i * 100
            ),
            start = start_preds,
            frequency = freq_x
          )

        }

        close(pb)
        snow::stopCluster(cl_SOCK)

      }
    }

    if (use_xreg) # with external regressors
    {
      if (!is.null(centers)) # with clustering
      {
        n_series_with_xreg_clusters <- n_series + n_xreg + centers
        preds_mean <- matrix(0, ncol = n_series_with_xreg_clusters, nrow = h)
        preds_upper <- matrix(0, ncol = n_series_with_xreg_clusters, nrow = h)
        preds_lower <- matrix(0, ncol = n_series_with_xreg_clusters, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      } else {
        n_series_with_xreg <- n_series + n_xreg
        preds_mean <- matrix(0, ncol = n_series_with_xreg, nrow = h)
        preds_upper <- matrix(0, ncol = n_series_with_xreg, nrow = h)
        preds_lower <- matrix(0, ncol = n_series_with_xreg, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      }
    } else { # without external regressors
      if (!is.null(centers))  # with clustering
      {
        n_series_with_clusters <- n_series + centers
        preds_mean <- matrix(0, ncol = n_series_with_clusters, nrow = h)
        preds_upper <- matrix(0, ncol = n_series_with_clusters, nrow = h)
        preds_lower <- matrix(0, ncol = n_series_with_clusters, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      } else { # without external regressors or clustering
        preds_mean <- matrix(0, ncol = n_series, nrow = h)
        preds_upper <- matrix(0, ncol = n_series, nrow = h)
        preds_lower <- matrix(0, ncol = n_series, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      }
    }

    type_aggregation <- match.arg(type_aggregation)

    for (j in 1:n_series)
    {
      sims_series_j <- sapply(1:B, function(i)
        sims[[i]][, j])
      preds_mean[, j] <- switch(type_aggregation,
                                mean = rowMeans(sims_series_j),
                                median = apply(sims_series_j, 1, median))
      preds_upper[, j] <-
        apply(sims_series_j, 1, function(x)
          quantile(x, probs = 1 - (1 - level / 100) / 2))
      preds_lower[, j] <-
        apply(sims_series_j, 1, function(x)
          quantile(x, probs = (1 - level / 100) / 2))
    }

    out <- list(
      mean = ts(
        data = preds_mean,
        start = start_preds,
        frequency = freq_x
      ),
      lower = ts(
        data = preds_lower,
        start = start_preds,
        frequency = freq_x
      ),
      upper = ts(
        data = preds_upper,
        start = start_preds,
        frequency = freq_x
      ),
      sims = sims,
      x = y,
      level = level,
      method = "ridge2",
      residuals = ts(fit_obj$resids,
                     start = start_x),
      loocv = fit_obj$loocv,
      weighted_loocv = fit_obj$weighted_loocv,
      loocv_per_series = fit_obj$loocv_per_series
    )

    if (use_xreg || use_clustering)
    {
      names_out <- names(out)
      for (i in 1:length(out))
      {
        try_delete_xreg <- try(delete_columns(out[[i]], "xreg_"),
                               silent = TRUE)
        if (!inherits(try_delete_xreg, "try-error") && !is.null(out[[i]]))
        {
          out[[i]] <- try_delete_xreg
        } else {
          if (identical(names_out[i], "sims")) # too much ifs man
          {
            # with simulations, it's a bit more tedious
            for (j in 1:B)
            {
              try_delete_xreg_sims <- try(delete_columns(out$sims[[j]], "xreg_"),
                                          silent = TRUE)
              if (!inherits(try_delete_xreg_sims, "try-error"))
              {
                out$sims[[j]] <- try_delete_xreg_sims
              }
            }
          }
        }
      }
    }

    if (!is.null(ym))
    {
      neutralized_sims <- lapply(series_names,
                                 function(series)
                                   neutralize(out, ym = ym,
                                              selected_series = series))
      names(neutralized_sims) <- series_names
      out$neutralized_sims <- neutralized_sims
    }

    if (dimensionality == "univariate")
    {
      names_out <- names(out)
      for (i in 1:length(out))
      {
        try_delete_trend <- try(delete_columns(out[[i]], "trend_univariate"),
                               silent = TRUE)
        if (!inherits(try_delete_trend, "try-error") && !is.null(out[[i]]))
        {
          out[[i]] <- try_delete_trend
        } else {
          if (identical(names_out[i], "sims")) # too much ifs man
          {
            # with simulations, it's a bit more tedious
            for (j in 1:B)
            {
              try_delete_trend_sims <- try(delete_columns(out$sims[[j]], "trend_univariate"),
                                          silent = TRUE)
              if (!inherits(try_delete_trend_sims, "try-error"))
              {
                out$sims[[j]] <- try_delete_trend_sims
              }
            }
          }
        }
      }
      return(structure(out, class = "forecast"))
    }

    return(structure(out, class = "mtsforecast"))
  }

  if (identical(type_pi, "rvinecopula"))
  {
    margins <- match.arg(margins)

    cl <- floor(min(max(cl, 0L), parallel::detectCores()))

    if (cl <= 1L)
    {
      # consider a loop with a progress bar
      sims <- lapply(1:B,
                     function(i)
                       ts(
                         fcast_ridge2_mts(
                           fit_obj,
                           h = h,
                           type_forecast = type_forecast,
                           type_simulation = type_pi,
                           margins = margins,
                           seed = seed + i * 100
                         ),
                         start = start_preds,
                         frequency = freq_x
                       ))
    } else {
      # consider a loop with a progress bar
      cluster <- parallel::makeCluster(getOption("cl.cores", cl))
      sims <-
        parallel::parLapply(
          cl = cluster,
          X = 1:B,
          fun = function(i)
            ts(
              fcast_ridge2_mts(
                fit_obj,
                h = h,
                type_forecast = type_forecast,
                type_simulation = type_pi,
                margins = margins,
                seed = seed + i * 100
              ),
              start = start_preds,
              frequency = freq_x
            )
        )
      parallel::stopCluster(cluster)
    }

    if (use_xreg) # with external regressors
    {
      if (!is.null(centers)) # with clustering
      {
        n_series_with_xreg_clusters <- n_series + n_xreg + centers
        preds_mean <- matrix(0, ncol = n_series_with_xreg_clusters, nrow = h)
        preds_upper <- matrix(0, ncol = n_series_with_xreg_clusters, nrow = h)
        preds_lower <- matrix(0, ncol = n_series_with_xreg_clusters, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      } else {
        n_series_with_xreg <- n_series + n_xreg
        preds_mean <- matrix(0, ncol = n_series_with_xreg, nrow = h)
        preds_upper <- matrix(0, ncol = n_series_with_xreg, nrow = h)
        preds_lower <- matrix(0, ncol = n_series_with_xreg, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      }
    } else { # without external regressors
      if (!is.null(centers))  # with clustering
      {
        n_series_with_clusters <- n_series + centers
        preds_mean <- matrix(0, ncol = n_series_with_clusters, nrow = h)
        preds_upper <- matrix(0, ncol = n_series_with_clusters, nrow = h)
        preds_lower <- matrix(0, ncol = n_series_with_clusters, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      } else { # without external regressors or clustering
        preds_mean <- matrix(0, ncol = n_series, nrow = h)
        preds_upper <- matrix(0, ncol = n_series, nrow = h)
        preds_lower <- matrix(0, ncol = n_series, nrow = h)
        colnames(preds_mean) <- series_names
        colnames(preds_upper) <- series_names
        colnames(preds_lower) <- series_names
      }
    }

    type_aggregation <- match.arg(type_aggregation)

    for (j in 1:n_series)
    {
      sims_series_j <- sapply(1:B, function(i)
        sims[[i]][, j])
      preds_mean[, j] <- switch(type_aggregation,
                                mean = rowMeans(sims_series_j),
                                median = apply(sims_series_j, 1, median))
      preds_upper[, j] <-
        apply(sims_series_j, 1, function(x)
          quantile(x, probs = 1 - (1 - level / 100) / 2))
      preds_lower[, j] <-
        apply(sims_series_j, 1, function(x)
          quantile(x, probs = (1 - level / 100) / 2))
    }

    out <- list(
      mean = ts(
        data = preds_mean,
        start = start_preds,
        frequency = freq_x
      ),
      lower = ts(
        data = preds_lower,
        start = start_preds,
        frequency = freq_x
      ),
      upper = ts(
        data = preds_upper,
        start = start_preds,
        frequency = freq_x
      ),
      sims = sims,
      x = y,
      level = level,
      method = "ridge2",
      residuals = ts(fit_obj$resids,
                     start = start_x),
      copula = fit_obj$params_distro,
      margins = margins,
      loocv = fit_obj$loocv,
      weighted_loocv = fit_obj$weighted_loocv,
      loocv_per_series = fit_obj$loocv_per_series
    )

    if (use_xreg || use_clustering) # refactor this
    {
      names_out <- names(out)
      for (i in 1:length(out))
      {
        try_delete_xreg <- try(delete_columns(out[[i]], "xreg_"),
                               silent = TRUE)
        if (!inherits(try_delete_xreg, "try-error") && !is.null(out[[i]]))
        {
          out[[i]] <- try_delete_xreg
        } else {
          if (identical(names_out[i], "sims")) # too much ifs man
          {
            # with simulations, it's a bit more tedious
            for (j in 1:B)
            {
              try_delete_xreg_sims <- try(delete_columns(out$sims[[j]], "xreg_"),
                                          silent = TRUE)
              if (!inherits(try_delete_xreg_sims, "try-error"))
              {
                out$sims[[j]] <- try_delete_xreg_sims
              }
            }
          }
        }
      }
    }

    if (!is.null(ym))
    {
      neutralized_sims <- lapply(series_names,
                                 function(series)
                                   neutralize(out, ym = ym,
                                              selected_series = series))
      names(neutralized_sims) <- series_names
      out$neutralized_sims <- neutralized_sims
    }

    if (dimensionality == "univariate") # refactor this
    {
      names_out <- names(out)
      for (i in 1:length(out))
      {
        try_delete_trend_univariate <- try(delete_columns(out[[i]], "trend_univariate"),
                               silent = TRUE)
        if (!inherits(try_delete_trend_univariate, "try-error") && !is.null(out[[i]]))
        {
          out[[i]] <- try_delete_trend_univariate
        } else {
          if (identical(names_out[i], "sims")) # too much ifs man
          {
            # with simulations, it's a bit more tedious
            for (j in 1:B)
            {
              try_delete_trend_univariate_sims <- try(delete_columns(out$sims[[j]], "trend_univariate"),
                                          silent = TRUE)
              if (!inherits(try_delete_trend_univariate_sims, "try-error"))
              {
                out$sims[[j]] <- try_delete_trend_univariate_sims
              }
            }
          }
        }
      }
      return(structure(out, class = "forecast"))
    }


    return(structure(out, class = "mtsforecast"))
  }

}
ridge2f <- compiler::cmpfun(ridge2f)


# Fitting function for ridge2
fit_ridge2_mts <- function(x,
                           lags = 1L,
                           nb_hidden = 5L,
                           nodes_sim = c("sobol", "halton", "unif"),
                           activ = c("relu", "sigmoid", "tanh",
                                     "leakyrelu", "elu", "linear"),
                           a = 0.01,
                           alpha = 0.5,
                           lambda_1 = 0.1,
                           lambda_2 = 0.1,
                           dropout = 0,
                           seed = 1,
                           type_pi = c(
                             "gaussian",
                             "bootstrap",
                             "blockbootstrap",
                             "movingblockbootstrap",
                             "rvinecopula"
                           ),
                           margins = c("gaussian",
                                       "empirical"),
                           hidden_layer_bias = FALSE,
                           ...)
{
  stopifnot(floor(nb_hidden) == nb_hidden)
  stopifnot(nb_hidden > 0)
  stopifnot(lambda_1 > 0 && lambda_2 > 0)

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_pi <- match.arg(type_pi)
  stopifnot(margins %in% c("gaussian", "empirical"))
  choice_margins <- match.arg(margins)
  margins <- switch(choice_margins, # this is convoluted
                    gaussian = "normal",
                    empirical = "empirical")

  series_names <- colnames(x)

  # /!\ because the ts object is in reverse order
  rev_x <- rev_matrix_cpp(x)

  # responses and regressors
  y_x <- create_train_inputs_cpp(rev_x, lags)

  # new predictors from original regressors
  list_regressors <- create_new_predictors(
    as.matrix(y_x$regressors),
    nb_hidden = nb_hidden,
    method = nodes_sim,
    activ = activ,
    a = a,
    seed = seed
  )
  new_regressors <- list_regressors$predictors

  # observed values, minus the lags (!)(beware)(!)
  observed_values <- y <- y_x$y

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(new_regressors)
  scaled_regressors <- x_scaled$res

  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  k_p <- lags * ncol(x)
  index <- 1:k_p

  # original predictors (scaled)
  X <- as.matrix(scaled_regressors[, index])

  # transformed predictors (scaled)
  Phi_X <- dropout_layer(scaled_regressors[,-index],
                         dropout = dropout, seed = seed)

  B <- crossprod(X) + lambda_1 * diag(k_p)
  C <- crossprod(Phi_X, X)
  D <- crossprod(Phi_X) + lambda_2 * diag(ncol(Phi_X))
  B_inv <- try(chol2inv(chol(B)), silent = TRUE)
  if (class(B_inv)[1] == "try-error") {
    B_inv <- solve(B) #my_ginv(B)
  }
  W <- C %*% B_inv
  S_mat <- D - tcrossprod(W, C)
  S_inv <- my_ginv(S_mat)
  Y <- S_inv %*% W
  inv <- rbind(cbind(B_inv + crossprod(W, Y), -t(Y)),
               cbind(-Y, S_inv))

  lscoef <- inv %*% crossprod(scaled_regressors, centered_y)
  colnames(lscoef) <- series_names
  rownames(lscoef) <- colnames(scaled_regressors)

  lsfit <- scaled_regressors %*% lscoef
  fitted_values <-
    rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                  ncol = ncol(lsfit)))
  resids <- rev_matrix_cpp(observed_values) - fitted_values
  colnames(resids) <- series_names

  params_distro <- NULL
  if (identical(type_pi, "rvinecopula")){
    params_distro <- try(select_residuals_dist(resids,
                                               uniformize = "ranks",
                                               distro = margins),
                         silent = TRUE)
    if (inherits(params_distro, "try-error"))
    {
      params_distro <- try(select_residuals_dist(resids,
                                                 uniformize = "ecdf",
                                                 distro = margins),
                           silent = TRUE)
    }
  }



  smoother_matrix <- scaled_regressors %*% tcrossprod(inv,
  scaled_regressors)
  loocv_matrix <- (resids / (1 - mean(diag(smoother_matrix))))^2
  loocv_per_series <- colMeans(loocv_matrix)
  names(loocv_per_series) <- series_names
  loocv_weights <- apply(resids, 2, sd)
  weighted_loocv <- mean(scale(loocv_matrix))
  loocv <- mean(loocv_matrix)

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      series_names = series_names,
      nb_hidden = nb_hidden,
      method = nodes_sim,
      w = list_regressors$w,
      activ = list_regressors$activ,
      activ_name = activ,
      a = a,
      lambda_1 = lambda_1,
      lambda_2 = lambda_2,
      seed = seed,
      nn_xm = list_regressors$xm,
      nn_xsd = list_regressors$xsd,
      ym = ym,
      xm = xm,
      scales = xsd,
      coef = lscoef,
      resids = resids,
      params_distro = params_distro,
      loocv = loocv,
      weighted_loocv = weighted_loocv,
      loocv_per_series = loocv_per_series
    )
  )
}


# Forecasting function for ridge2
fcast_ridge2_mts <- function(fit_obj,
                             h = 5,
                             type_forecast = c("recursive", "direct"),
                             level = 95,
                             type_simulation = c("none", "rvinecopula"), # other than bootstrap
                             margins = c("gaussian", "empirical"),
                             bootstrap = FALSE,
                             type_bootstrap = c("bootstrap",
                                                "blockbootstrap",
                                                "movingblockbootstrap"),
                             block_length = NULL,
                             seed = 123)
{
  type_forecast <- match.arg(type_forecast)
  type_simulation <- match.arg(type_simulation)
  margins <- switch(match.arg(margins),
                   gaussian = "normal",
                   empirical = "empirical")

  if (identical(bootstrap, FALSE))
  {

    if(identical(type_simulation, "none")) # other than bootstrap
    {

      # 1 - recursive forecasts (bootstrap == FALSE) -------------------------------------------------

      # recursive forecasts
      if (identical(type_forecast, "recursive"))
      {
        # observed values (minus lagged) in decreasing order (most recent first)
        y_mat <- rbind(matrix(NA, nrow = h, ncol = ncol(fit_obj$y)),
                       fit_obj$y)
        y <- fit_obj$y
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          newx <- cbind(newx, matrix(g(
            my_scale(newx, xm = nn_xm,
                     xsd = nn_xsd) %*% w
          ), nrow = 1))
          y_mat[h - i + 1, ] <- predict_myridge(fit_obj, newx = newx)
          y <- y_mat[complete.cases(y_mat), ]
        }
      }

      # 2 - direct forecasts (bootstrap == FALSE) -------------------------------------------------

      # direct forecasts
      if (identical(type_forecast, "direct"))
      {
        # observed values (minus lagged)
        y <- fit_obj$y
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          newx <- cbind(newx, matrix(g(
            my_scale(newx, xm = nn_xm,
                     xsd = nn_xsd) %*% w
          ), nrow = 1))


          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <-
            rbind(fit_obj$x, preds)[-1,] # same window length as x
          fit_obj <-
            fit_ridge2_mts(
              x = newtrainingx,
              lags = fit_obj$lags,
              nb_hidden = fit_obj$nb_hidden,
              nodes_sim = fit_obj$method,
              activ = fit_obj$activ_name,
              a = fit_obj$a,
              lambda_1 = fit_obj$lambda_1,
              lambda_2 = fit_obj$lambda_2,
              seed = fit_obj$seed
            )
        }
      }

      res2 <- rev_matrix_cpp(y)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n, ]
      colnames(res) <- fit_obj$series_names
      return(res)

    }

    if(identical(type_simulation, "rvinecopula"))
    {
      # 1 - recursive forecasts (bootstrap == FALSE, type_simulation == "rvinecopula") -------------------------------------------------

      # recursive forecasts
      if (identical(type_forecast, "recursive"))
      {
        # observed values (minus lagged) in decreasing order (most recent first)
        y <- fit_obj$y
        y_mat <- rbind(matrix(NA, nrow = h, ncol = ncol(fit_obj$y)),
                       fit_obj$y)
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          newx <- cbind(newx, matrix(g(
            my_scale(newx, xm = nn_xm,
                     xsd = nn_xsd) %*% w
          ), nrow = 1))
          y_mat[h - i + 1, ] <- predict_myridge(fit_obj,
                                                newx = newx)
          y <- y_mat[complete.cases(y_mat), ]
        }
      }

      # 2 - direct forecasts (bootstrap == FALSE, type_simulation == "rvinecopula") -------------------------------------------------
      # direct forecasts
      if (identical(type_forecast, "direct"))
      {
        # observed values (minus lagged)
        y <- fit_obj$y
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          newx <- cbind(newx, matrix(g(
            my_scale(newx, xm = nn_xm,
                     xsd = nn_xsd) %*% w
          ), nrow = 1))


          preds <- predict_myridge(fit_obj,
                                   newx = newx)

          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <-
            rbind(fit_obj$x, preds)[-1,] # same window length as x
          fit_obj <-
            fit_ridge2_mts(
              x = newtrainingx,
              lags = fit_obj$lags,
              nb_hidden = fit_obj$nb_hidden,
              nodes_sim = fit_obj$method,
              activ = fit_obj$activ_name,
              a = fit_obj$a,
              lambda_1 = fit_obj$lambda_1,
              lambda_2 = fit_obj$lambda_2,
              seed = fit_obj$seed
            )
        }
      }

      residuals_simulations <- simulate_rvine(fit_obj,
                                              RVM_U = fit_obj$params_distro$RVM_U,
                                              h = h,
                                              seed = seed,
                                              tests = FALSE)

      res2 <- rev_matrix_cpp(y)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n, ] + as.matrix(residuals_simulations)
      colnames(res) <- fit_obj$series_names
      return(res)
    }

  } else { # if (bootstrap == TRUE)
    # if bootstrap == TRUE
    type_bootstrap <- match.arg(type_bootstrap)

    # observed values (minus lagged) in decreasing order (most recent first)
    y <- fit_obj$y
    freq_y <- frequency(y)
    if (nrow(y) <= 2 * freq_y)
      freq_y <- 1L

    if (type_bootstrap %in% c("blockbootstrap", "movingblockbootstrap"))
    {
      if (is.null(block_length)) {
        block_length <- ifelse(freq_y > 1, 2 * freq_y, min(8, floor(nrow(y) / 2)))
      }
    }

    # sampling from the residuals independently or in blocks
    set.seed(seed)

    if (type_bootstrap %in% c("bootstrap", "blockbootstrap", "movingblockbootstrap"))
    {
      idx <- switch(
        type_bootstrap,
        bootstrap = sample.int(
          n = nrow(fit_obj$resids),
          size = h,
          replace = TRUE
        ),
        blockbootstrap = mbb(
          r = fit_obj$resids,
          n = h,
          b = block_length,
          return_indices = TRUE,
          seed = seed
        ), # in utils.R
        movingblockbootstrap = mbb2(
          r = fit_obj$resids,
          n = h,
          b = block_length,
          return_indices = TRUE,
          seed = seed
        ) # in utils.R
      )
    }


    # 1 - recursive forecasts (bootstrap == TRUE) -------------------------------------------------

    # recursive forecasts
    if (identical(type_forecast, "recursive"))
    {
      y_mat <- rbind(matrix(NA, nrow = h, ncol = ncol(fit_obj$y)),
                     fit_obj$y)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd

      w <- fit_obj$w
      g <- fit_obj$activ

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)

        newx <- cbind(newx, matrix(g(
          my_scale(newx, xm = nn_xm,
                   xsd = nn_xsd) %*% w
        ), nrow = 1))
        y_mat[h - i + 1, ] <-
          predict_myridge(fit_obj, newx = newx) + fit_obj$resids[idx[i],]
        y <- y_mat[complete.cases(y_mat), ]
      }
    }

    # 2 - direct forecasts (bootstrap == TRUE) -------------------------------------------------

    # direct forecasts
    if (identical(type_forecast, "direct"))
    {
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd

      w <- fit_obj$w
      g <- fit_obj$activ

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)

        newx <- cbind(newx, matrix(g(
          my_scale(newx, xm = nn_xm,
                   xsd = nn_xsd) %*% w
        ), nrow = 1))


        preds <- predict_myridge(fit_obj, newx = newx) + fit_obj$resids[idx[i],]
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1,] # same window length as x
        fit_obj <- fit_ridge2_mts(x = newtrainingx,
                                  lags = fit_obj$lags,
                                  nb_hidden = fit_obj$nb_hidden,
                                  nodes_sim = fit_obj$method,
                                  activ = fit_obj$activ_name,
                                  a = fit_obj$a,
                                  lambda_1 = fit_obj$lambda_1,
                                  lambda_2 = fit_obj$lambda_2,
                                  seed = fit_obj$seed)
      }
    }

    res2 <- rev_matrix_cpp(y)
    n <- nrow(res2)
    res <- res2[(n - h + 1):n, ]
    colnames(res) <- fit_obj$series_names
    return(res)
  }

}
