



#' Ridge2 model
#'
#' Random Vector functional link network model with 2 regularization parameters
#'
#' @param y A multivariate time series of class \code{ts} (preferred) or a \code{matrix}
#' @param xreg External regressors. A data.frame (preferred) or a \code{matrix}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags
#' @param nb_hidden Number of nodes in hidden layer
#' @param nodes_sim Type of simulation for nodes in the hidden layer
#' @param activ Activation function
#' @param a hyperparameter for activation function "leakyrelu", "elu"
#' @param lambda_1 Regularization parameter for original predictors
#' @param lambda_2 Regularization parameter for transformed predictors
#' @param dropout dropout regularization parameter (dropping nodes in hidden layer)
#' @param seed Reproducibility seed for `nodes_sim == unif`
#' @param type_forecast Recursive or direct forecast
#' @param type_pi type of prediction interval currently "gaussian" or (independent) "bootstrap" or (circular) "blockbootstrap"
#' @param block_length length of block for circular block bootstrap (\code{type_pi == 'blockbootstrap'})
#' @param seed reproducibility seed for \code{type_pi == 'bootstrap'} or \code{type_pi == 'blockbootstrap'}
#' @param B Number of bootstrap replications for \code{type_pi == 'bootstrap'} or \code{type_pi == 'blockbootstrap'}
#' @param cl an integer; the number of clusters for parallel execution, for \code{type_pi == 'bootstrap'} or \code{type_pi == 'blockbootstrap'}
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
ridge2f <- function(y,
                    xreg = NULL,
                    h = 5,
                    level = 95,
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
                    type_pi = c("gaussian", "bootstrap", "blockbootstrap"),
                    block_length = NULL,
                    seed = 1,
                    B = 100L,
                    cl = 1L)
{
  stopifnot(!is.null(ncol(y)))

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

  if (!is.null(xreg))
  {
    if (is.null(ncol(xreg)))
      xreg <- as.matrix(xreg)

    stopifnot(nrow(xreg) == nrow(y))

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

  series_names <- colnames(y)
  type_pi <- match.arg(type_pi)
  freq_x <- frequency(y)
  start_x <- start(y)
  start_preds <- tsp(y)[2] + 1 / freq_x

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_forecast <- match.arg(type_forecast)

  # Fitting a regularized regression  model to multiple time series
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
    seed = seed
  )

  if (type_pi == "gaussian")
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
      residuals = fit_obj$resids,
      coefficients = fit_obj$coef
    )

    if (use_xreg)
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

    return(structure(out, class = "mtsforecast"))
  }

  if (type_pi %in% c("bootstrap", "blockbootstrap"))
  {
    if (cl <= 1L)
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
    } else {
      stopifnot(is.numeric(cl))
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
    }

    if (use_xreg)
    {
      n_series_with_xreg <- n_series + n_xreg
      preds_mean <- matrix(0, ncol = n_series_with_xreg, nrow = h)
      preds_upper <- matrix(0, ncol = n_series_with_xreg, nrow = h)
      preds_lower <- matrix(0, ncol = n_series_with_xreg, nrow = h)
      colnames(preds_mean) <- series_names
      colnames(preds_upper) <- series_names
      colnames(preds_lower) <- series_names
    } else {
      preds_mean <- matrix(0, ncol = n_series, nrow = h)
      preds_upper <- matrix(0, ncol = n_series, nrow = h)
      preds_lower <- matrix(0, ncol = n_series, nrow = h)
      colnames(preds_mean) <- series_names
      colnames(preds_upper) <- series_names
      colnames(preds_lower) <- series_names
    }

    for (j in 1:n_series)
    {
      sims_series_j <- sapply(1:B, function(i)
        sims[[i]][, j])
      preds_mean[, j] <- rowMeans(sims_series_j)
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
      residuals = fit_obj$resids
    )

    if (use_xreg)
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
                           seed = 1)
{
  stopifnot(floor(nb_hidden) == nb_hidden)
  stopifnot(nb_hidden > 0)
  stopifnot(lambda_1 > 0 && lambda_2 > 0)

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

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
      resids = resids
    )
  )
}



# Forecasting function for ridge2
fcast_ridge2_mts <- function(fit_obj,
                             h = 5,
                             type_forecast = c("recursive", "direct"),
                             level = 95,
                             bootstrap = FALSE,
                             type_bootstrap = c("bootstrap", "blockbootstrap"),
                             block_length = NULL,
                             seed = 123)
{
  type_forecast <- match.arg(type_forecast)

  if (bootstrap == FALSE)
  {
    # 1 - recursive forecasts (bootstrap == FALSE) -------------------------------------------------

    # recursive forecasts
    if (type_forecast == "recursive")
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
    if (type_forecast == "direct")
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

  } else {

    # if bootstrap == TRUE
    type_bootstrap <- match.arg(type_bootstrap)

    # sampling from the residuals independently or in blocks
    if (type_bootstrap %in% c("bootstrap", "blockbootstrap"))
    {
      if (type_bootstrap == "blockbootstrap")
        stopifnot(!is.null(block_length)) # use rlang::abort

      set.seed(seed)
      idx <- switch(type_bootstrap,
                    bootstrap = sample.int(n = nrow(fit_obj$resids),
                                           size = h,
                                           replace = TRUE),
                    blockbootstrap = mbb(
                      r = fit_obj$resids,
                      n = h,
                      b = block_length,
                      return_indices = TRUE,
                      seed = seed
                    ))
      # cat("idx: ", idx, "\n")
    }

    # 1 - recursive forecasts (bootstrap == TRUE) -------------------------------------------------

    # recursive forecasts
    if (type_forecast == "recursive")
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
        y_mat[h - i + 1, ] <-
          predict_myridge(fit_obj, newx = newx) + fit_obj$resids[idx[i],]
        y <- y_mat[complete.cases(y_mat), ]
      }
    }

    # 2 - direct forecasts (bootstrap == TRUE) -------------------------------------------------

    # direct forecasts
    if (type_forecast == "direct")
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


        preds <-
          predict_myridge(fit_obj, newx = newx) + fit_obj$resids[idx[i],]
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

}
