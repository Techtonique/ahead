
#' Ridge2 model
#'
#' Random Vector functional link network model with 2 regularization parameters
#'
#' @param y A multivariate time series of class \code{ts}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags
#' @param nb_hidden Number of nodes in hidden layer
#' @param nodes_sim Type of simulation for nodes in the hidden layer
#' @param activ Activation function
#' @param a hyperparameter for activation function "leakyrelu", "elu"
#' @param lambda_1 Regularization parameter for original predictors
#' @param lambda_2 Regularization parameter for transformed predictors
#' @param seed Reproducibility seed for `nodes_sim == unif`
#' @param type_forecast Recursive or direct forecast
#' @param type_pi currently "gaussian" or "bootstrap"
#' @param seed reproducibility seed for \code{type_pi == 'bootstrap'}
#' @param B Number of bootstrap replications for \code{type_pi == 'bootstrap'}
#'
#' @return An object of class "mtsforecast"; a list containing the following elements:
#'
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{x}{The original time series}
#' \item{residuals}{Residuals from the fitted model}
#' \item{sims}{Model simulations for \code{type_pi == bootstrap}}
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
#' res <- ahead::ridge2f(fpp::usconsumption, h=20, lags=2,
#' lambda_2=1)
#' par(mfrow=c(1, 2))
#' plot(res, "income")
#' plot(res, "consumption")
#'
ridge2f <- function(y, h = 5,
                    level = 95,
                    lags = 1,
                    nb_hidden = 5,
                    nodes_sim = c("sobol", "halton", "unif"),
                    activ = c("relu", "sigmoid", "tanh",
                              "leakyrelu", "elu", "linear"),
                    a = 0.01,
                    lambda_1 = 0.1,
                    lambda_2 = 0.1,
                    type_forecast = c("recursive", "direct"),
                    type_pi = c("gaussian", "bootstrap"),
                    seed = 1, B = 100)
{
  if (!is.ts(y))
  {
    y <- ts(y)
  }
  stopifnot(!is.null(ncol(y)))
  n_series <- ncol(y)
  series_names <- colnames(y)
  type_pi <- match.arg(type_pi)
  freq_x <- frequency(y)
  start_x <- start(y)
  start_preds <- tsp(y)[2] + 1/freq_x

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
    seed = seed
  )

  if (type_pi == "gaussian")
  {
    preds <- ts(data = fcast_ridge2_mts(
      fit_obj,
      h = h,
      type_forecast = type_forecast,
      level = level),
      start = start_preds, frequency = freq_x)

    # strong gaussian assumption
    sds <- apply(fit_obj$resids, 2, sd)
    qtsd <- -qnorm(0.5 - level/200)*sds
    qtsd <- t(replicate(nrow(preds), qtsd))

    # Forecast from fit_obj
    out <- list(mean = preds,
                lower = preds - qtsd,
                upper = preds + qtsd,
                sims = NULL,
                x = y,
                level = level,
                method = "ridge2",
                residuals = fit_obj$resids)

    # cat("Point Forecast", "\n")
    # print(out$mean)
    # cat("\n")
    #
    # cat("Lo", level, "\n")
    # print(out$lower)
    # cat("\n")
    #
    # cat("Hi", level, "\n")
    # print(out$upper)
    # cat("\n")

    return(structure(out, class = "mtsforecast"))
  }

  if (type_pi == "bootstrap")
  {
    sims <- lapply(1:B,
                    function(i) ts(fcast_ridge2_mts(fit_obj,
                    h = h,
                    type_forecast = type_forecast,
                    bootstrap = TRUE, seed = seed + i*100),
                    start = start_preds,
                    frequency = freq_x))

    `%op%` <- foreach::`%do%`

    j <- NULL

    preds_mean <- foreach::foreach(j = 1:n_series,
                                  .combine = cbind)%op%{
                                  rowMeans(sapply(1:B, function(i) sims[[i]][, j]))
                                  }
    colnames(preds_mean) <- series_names

    preds_upper <- foreach::foreach(j = 1:n_series,
                                   .combine = cbind)%op%{
                                    preds_j_upper <- sapply(1:B, function(i) sims[[i]][, j])
                                    apply(preds_j_upper, 1, function(x) quantile(x, probs = level/100))
                                   }
    colnames(preds_upper) <- series_names

    preds_lower <- foreach::foreach(j = 1:n_series,
                                    .combine = cbind)%op%{
                                     preds_j_lower <- sapply(1:B, function(i) sims[[i]][, j])
                                     apply(preds_j_lower, 1, function(x) quantile(x, probs = 1 - level/100))
                                    }
    colnames(preds_lower) <- series_names

    # colnames
    out <- list(mean = ts(data = preds_mean,
                          start = start_preds, frequency = freq_x),
                lower = ts(data = preds_lower,
                           start = start_preds, frequency = freq_x),
                upper = ts(data = preds_upper,
                           start = start_preds, frequency = freq_x),
                sims = sims,
                x = y,
                level = level,
                method = "ridge2",
                residuals = fit_obj$resids)

    # cat("Point Forecast", "\n")
    # print(out$mean)
    # cat("\n")
    #
    # cat("Lo", level, "\n")
    # print(out$lower)
    # cat("\n")
    #
    # cat("Hi", level, "\n")
    # print(out$upper)
    # cat("\n")

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
                           seed = 1)
{
  stopifnot(floor(nb_hidden)==nb_hidden)
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
  Phi_X <- scaled_regressors[, -index]

  B <- crossprod(X) + lambda_1 * diag(k_p)
  C <- crossprod(Phi_X, X)
  D <- crossprod(Phi_X) + lambda_2 * diag(ncol(Phi_X))
  B_inv <- solve(B) #my_ginv(B)
  W <- C %*% B_inv
  S_mat <- D - tcrossprod(W, C)
  S_inv <- my_ginv(S_mat)
  Y <- S_inv %*% W
  inv <- rbind(cbind(B_inv + crossprod(W, Y),-t(Y)),
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
                          level = 95, # for gaussian p.i
                          bootstrap = FALSE, # for bootstrap p.i
                          seed=123) # for bootstrap p.i
{
  type_forecast <- match.arg(type_forecast)

  if(bootstrap == FALSE)
  {
    # 1 - recursive forecasts -------------------------------------------------

    # recursive forecasts
    if (type_forecast == "recursive")
    {
      # observed values (minus lagged) in decreasing order (most recent first)
      y <- fit_obj$y
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd

      w <- fit_obj$w
      g <- fit_obj$activ

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)

        newx <- cbind(newx, matrix(g(my_scale(
          newx, xm = nn_xm,
          xsd = nn_xsd
        ) %*% w), nrow = 1))

        preds <- predict_myridge(fit_obj, newx = newx)
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    # 2 - direct forecasts -------------------------------------------------

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

        newx <- cbind(newx, matrix(g(my_scale(
          newx, xm = nn_xm,
          xsd = nn_xsd
        ) %*% w), nrow = 1))


        preds <- predict_myridge(fit_obj, newx = newx)
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
        fit_obj <- fit_ridge2_mts(x = newtrainingx, lags = fit_obj$lags,
                                  nb_hidden = fit_obj$nb_hidden,
                                  nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                  a = fit_obj$a, lambda_1 = fit_obj$lambda_1, lambda_2 = fit_obj$lambda_2,
                                  seed = fit_obj$seed)
      }
    }

    res2 <- rev_matrix_cpp(y)
    n <- nrow(res2)
    res <- res2[(n - h + 1):n,]
    colnames(res) <- fit_obj$series_names
    return(res)

  } else { # if bootstrap == TRUE

    # sampling from the residuals
    set.seed(seed)
    idx <- sample.int(n=nrow(fit_obj$resids), size=h, replace = TRUE)

    # 1 - recursive forecasts -------------------------------------------------

    # recursive forecasts
    if (type_forecast == "recursive")
    {
      # observed values (minus lagged) in decreasing order (most recent first)
      y <- fit_obj$y
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd

      w <- fit_obj$w
      g <- fit_obj$activ

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)

        newx <- cbind(newx, matrix(g(my_scale(
          newx, xm = nn_xm,
          xsd = nn_xsd
        ) %*% w), nrow = 1))

        preds <- predict_myridge(fit_obj, newx = newx) + fit_obj$resids[idx[i], ]
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    # 2 - direct forecasts -------------------------------------------------

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

        newx <- cbind(newx, matrix(g(my_scale(
          newx, xm = nn_xm,
          xsd = nn_xsd
        ) %*% w), nrow = 1))


        preds <- predict_myridge(fit_obj, newx = newx) + fit_obj$resids[idx[i], ]
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
        fit_obj <- fit_ridge2_mts(x = newtrainingx, lags = fit_obj$lags,
                                  nb_hidden = fit_obj$nb_hidden,
                                  nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                  a = fit_obj$a, lambda_1 = fit_obj$lambda_1, lambda_2 = fit_obj$lambda_2,
                                  seed = fit_obj$seed)
      }
    }

    res2 <- rev_matrix_cpp(y)
    n <- nrow(res2)
    res <- res2[(n - h + 1):n,]
    colnames(res) <- fit_obj$series_names
    return(res)
  }

}
