# In alphabetical order

# create new predictors -----
create_new_predictors <- function(x,
                                  nb_hidden = 5,
                                  hidden_layer_bias = FALSE,
                                  method = c("sobol", "halton", "unif"),
                                  activ = c("relu", "sigmoid", "tanh",
                                            "leakyrelu", "elu", "linear"),
                                  a = 0.01,
                                  seed = 123)
{
  n <- nrow(x)

  if (nb_hidden > 0)
  {
    p <- ncol(x)
    method <- match.arg(method)

    # Activation function
    g <- switch(
      match.arg(activ),
      "relu" = function(x)
        x * (x > 0),
      "sigmoid" = function(x)
        1 / (1 + exp(-x)),
      "tanh" = function(x)
        tanh(x),
      "leakyrelu" = function(x)
        x * (x > 0) + a * x * (x <= 0),
      "elu" = function(x)
        x * (x >= 0) + a * (exp(x) - 1) * (x < 0),
      "linear" = function(x)
        x
    )

    if (hidden_layer_bias == FALSE)
    {
      # used for columns sample and for 'method == unif'
      set.seed(seed + 1)
      w <- remove_zero_cols(switch(
        method,
        "sobol" = 2 * t(randtoolbox::sobol(nb_hidden + 1, p)) - 1,
        "halton" = 2 * t(randtoolbox::halton(nb_hidden, p)) - 1,
        "unif" = matrix(
          runif(nb_hidden * p, min = -1, max = 1),
          nrow = p,
          ncol = nb_hidden
        )
      ))
      scaled_x <- my_scale(x)
      hidden_layer_obj <- remove_zero_cols(g(scaled_x$res %*% w),
                                           with_index = TRUE)
      hidden_layer <- hidden_layer_obj$mat

    } else {
      # hidden_layer_bias == TRUE
      pp <- p + 1
      # used for columns sample and for 'method == unif'
      set.seed(seed + 1)
      w <- remove_zero_cols(switch(
        method,
        "sobol" = 2 * t(randtoolbox::sobol(nb_hidden + 1, pp)) - 1,
        "halton" = 2 * t(randtoolbox::halton(nb_hidden, pp)) - 1,
        "unif" = matrix(
          runif(nb_hidden * pp, min = -1, max = 1),
          nrow = pp,
          ncol = nb_hidden
        )
      ))

      scaled_x <- my_scale(x)
      hidden_layer_obj <-
        remove_zero_cols(g(cbind(1, scaled_x$res) %*% w),
                         with_index = TRUE)
      hidden_layer <- hidden_layer_obj$mat
    }

    res <- cbind(x, hidden_layer)
    nb_nodes <- ncol(hidden_layer)
    if (!is.null(nb_nodes))
      colnames(res) <-
      c(paste0("x", 1:p), # maybe use the real names
        paste0("h", 1:nb_nodes))


    # if nb_hidden > 0 && (nb_predictors >= 2 && col_sample < 1)
    return(
      list(
        activ = g,
        xm = scaled_x$xm,
        xsd = scaled_x$xsd,
        w = w,
        predictors = res,
        hidden_layer_index = hidden_layer_obj$index
      )
    )
  } else {
    # if nb_hidden <= 0
    scaled_x <- my_scale(x)
    return(
      list(
        xm = scaled_x$xm,
        xsd = scaled_x$xsd,
        predictors = x,
        hidden_layer_index = hidden_layer_obj$index
      )
    )
  }
}

# create trend and seasonality features -----

#' Create trend and seasonality features for univariate time series
#'
#' @param y a univariate time series
#'
#' @return a vector or matrix of features
#' @export
#'
#' @examples
#'
#' y <- ts(rnorm(100), start = 1, frequency = 12)
#' createtrendseason(y)
#'
#' createtrendseason(USAccDeaths)
#'
createtrendseason <- function(y) {
  seq_along_y <- seq_along(y)
  frequency_y <- frequency(y)
  theta <- 2 * pi * seq_along_y/frequency_y
  sin_season <- sin(theta)
  cos_season <- cos(theta)

  if(all(cos_season == 1) || sum(sin_season^2) < .Machine$double.eps) {
    return(ts(seq_along(y),
              start = start(y)))
  }
  result_data <- data.frame(
    sin_season = sin_season,
    cos_season = cos_season
  )
  res <- ts(as.matrix(cbind(1:length(y),
                            result_data)),
            start = start(y),
            frequency = frequency(y))
  colnames(res) <- c("trend", "sin_season", "cos_season")
  return(res)
}

# delete columns using a string pattern -----
delete_columns <- function(x, pattern)
{
  x[, !grepl(pattern = pattern, x = colnames(x))]
}

# dropout regularization -----
dropout_layer <- function(X, dropout = 0, seed = 123)
{
  stopifnot(dropout <= 0.8)
  stopifnot(dropout >= 0)
  if (dropout == 0)
  {
    return(X)
  } else {
    n_rows <- dim(X)[1]
    n_columns <- dim(X)[2]
    set.seed(seed)
    mask <- (matrix(
      runif(n_rows * n_columns),
      nrow = n_rows,
      ncol = n_columns
    ) > dropout)
    return(X * mask / (1 - dropout))
  }
}

# MASS::fitdistr, but with stats::nlminb -----
#'
#' @export
#'
fitdistr_ahead <- function (x, densfun, start = NULL, seed = 123, ...)
{
  myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))
  mydt <- function(x, m, s, df, log) stats::dt((x - m)/s, df, log = TRUE) - log(s)
  Call <- match.call(expand.dots = TRUE)
  if (missing(start))
    start <- NULL
  dots <- names(list(...))
  dots <- dots[!is.element(dots, c("upper", "lower"))]
  if (missing(x) || length(x) == 0L || mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x)))
    stop("'x' contains missing or infinite values")
  if (missing(densfun) || !(is.function(densfun) || is.character(densfun)))
    stop("'densfun' must be supplied as a function or name")
  control <- list()
  n <- length(x)
  if (is.character(densfun)) {
    distname <- tolower(densfun)
    densfun <- switch(distname, beta = stats::dbeta, cauchy = stats::dcauchy,
                      `chi-squared` = stats::dchisq, exponential = stats::dexp, f = stats::df,
                      gamma = stats::dgamma, geometric = stats::dgeom, `log-normal` = stats::dlnorm,
                      lognormal = stats::dlnorm, logistic = stats::dlogis, `negative binomial` = stats::dnbinom,
                      normal = stats::dnorm, poisson = stats::dpois, t = stats::dt,
                      weibull = stats::dweibull, NULL)
    if (is.null(densfun))
      stop("unsupported distribution")
    if (distname %in% c("lognormal", "log-normal")) {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "log-Normal"), domain = NA)
      if (any(x <= 0))
        stop("need positive values to fit a log-Normal")
      lx <- log(x)
      sd0 <- sqrt((n - 1)/n) * sd(lx)
      mx <- mean(lx)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
                   dimnames = list(names(sds), names(sds)))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(stats::dlnorm(x, mx,
                                                                  sd0, log = TRUE))), class = "fitdistr"))
    }
    if (distname == "normal") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "Normal"), domain = NA)
      sd0 <- sqrt((n - 1)/n) * sd(x)
      mx <- mean(x)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("mean", "sd")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
                   dimnames = list(names(sds), names(sds)))
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(stats::dnorm(x, mx, sd0,
                                                                 log = TRUE))), class = "fitdistr"))
    }
    if (distname == "poisson") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "Poisson"), domain = NA)
      estimate <- mean(x)
      sds <- sqrt(estimate/n)
      names(estimate) <- names(sds) <- "lambda"
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("lambda",
                                                              "lambda"))
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(stats::dpois(x, estimate,
                                                                 log = TRUE))), class = "fitdistr"))
    }
    if (distname == "exponential") {
      if (any(x < 0))
        stop("Exponential values must be >= 0")
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "exponential"), domain = NA)
      estimate <- 1/mean(x)
      sds <- estimate/sqrt(n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("rate",
                                                              "rate"))
      names(estimate) <- names(sds) <- "rate"
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(stats::dexp(x, estimate,
                                                                log = TRUE))), class = "fitdistr"))
    }
    if (distname == "geometric") {
      if (!is.null(start))
        stop(gettextf("supplying pars for the %s distribution is not supported",
                      "geometric"), domain = NA)
      estimate <- 1/(1 + mean(x))
      sds <- estimate * sqrt((1 - estimate)/n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("prob",
                                                              "prob"))
      names(estimate) <- names(sds) <- "prob"
      return(structure(list(estimate = estimate, sd = sds,
                            vcov = vc, n = n, loglik = sum(stats::dgeom(x, estimate,
                                                                 log = TRUE))), class = "fitdistr"))
    }
    if (distname == "weibull" && is.null(start)) {
      if (any(x <= 0))
        stop("Weibull values must be > 0")
      lx <- log(x)
      m <- mean(lx)
      v <- stats::var(lx)
      shape <- 1.2/sqrt(v)
      scale <- exp(m + 0.572/shape)
      start <- list(shape = shape, scale = scale)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "gamma" && is.null(start)) {
      if (any(x < 0))
        stop("gamma values must be >= 0")
      m <- mean(x)
      v <- stats::var(x)
      start <- list(shape = m^2/v, rate = m/v)
      start <- start[!is.element(names(start), dots)]
      control <- list(parscale = c(1, start$rate))
    }
    if (distname == "negative binomial" && is.null(start)) {
      m <- mean(x)
      v <- stats::var(x)
      size <- if (v > m)
        m^2/(v - m)
      else 100
      start <- list(size = size, mu = m)
      start <- start[!is.element(names(start), dots)]
    }
    if (is.element(distname, c("cauchy", "logistic")) &&
        is.null(start)) {
      start <- list(location = median(x), scale = stats::IQR(x)/2)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "t" && is.null(start)) {
      start <- list(m = median(x), s = stats::IQR(x)/2, df = 10)
      start <- start[!is.element(names(start), dots)]
    }
  }
  if (is.null(start) || !is.list(start))
    stop("'start' must be a named list")
  nm <- names(start)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m)))
    stop("'start' specifies names which are not arguments to 'densfun'")
  formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
  dens <- function(parm, x, ...) densfun(x, parm, ...)
  if ((l <- length(nm)) > 1L)
    body(dens) <- parse(text = paste("densfun(x,", paste("parm[",
                                                         1L:l, "]", collapse = ", "), ", ...)"))

  Call[[1L]] <- quote(stats::nlminb)
  Call$densfun <- Call$start <- NULL
  Call$x <- x
  Call$start <- start

  Call$objective <- if ("log" %in% args)
    mylogfn
  else myfn

  Call$hessian <- TRUE # in MASS::fitdistr
  # Call$hessian <- NULL
  if (length(control))
    Call$control <- control

  res <- eval.parent(Call)

  return(list(estimate = res$par,
              convergence = res$convergence,
              objective = res$objective))
}

# clustering matrix -----
get_clusters <- function(x,
                         centers,
                         type_clustering = c("kmeans", "hclust"),
                         start = NULL,
                         frequency = NULL,
                         seed = 123,
                         ...)
{
  stopifnot(!missing(x)) # use rlang::abort

  stopifnot(!missing(centers)) # use rlang::abort

  # /!\ important
  x_scaled <- scale(x = x,
                    scale = TRUE,
                    center = TRUE)[,]

  type_clustering <- match.arg(type_clustering)

  set.seed(seed)

  df_clusters <- switch(
    type_clustering,
    kmeans = data.frame(stats::kmeans(x_scaled,
                                      centers = centers, ...)$cluster),
    hclust = data.frame(stats::cutree(stats::hclust(
      stats::dist(x_scaled, ...), ...
    ),
    k = centers))
  )

  df_clusters[[1]] <- as.factor(df_clusters[[1]])

  matrix_clusters <- stats::model.matrix( ~ -1 + ., df_clusters)

  rownames(matrix_clusters) <- NULL

  colnames(matrix_clusters) <- NULL

  matrix_clusters <- matrix_clusters[, ]

  if (!is.null(start) && !is.null(frequency))
  {
    cluster_ts <- ts(matrix_clusters,
                     start = start, frequency = frequency)

    colnames(cluster_ts) <-
      paste0("xreg_cluster", 1:ncol(cluster_ts))

    return(cluster_ts)
  } else {
    colnames(matrix_clusters) <- paste0("xreg_cluster", 1:ncol(matrix_clusters))
    return(matrix_clusters)
  }
}

# forecasting from a fitted object -----
fcast_obj_mts <- function(fit_obj,
                          h = 5,
                          type_pi = "none",
                          type_forecast = c("recursive", "direct"),
                          level = 95)
{
  # Fit method
  fit_method <- fit_obj$class_res
  # Method for forecasting the residuals
  stopifnot(!is.na(match(
    type_pi, c("none", "arima",
               "ets", "theta",
               "VAR")
  )))
  type_forecast <- match.arg(type_forecast)
  # index of chosen predictors (TRUE if all, NumericVector if col_sample < 1)
  col_index <- fit_obj$col_index
  # non-zero columns index in hidden layer
  hidden_layer_index <- fit_obj$hidden_layer_index

  # 1 - recursive forecasts -------------------------------------------------

  # recursive forecasts
  if (type_forecast == "recursive")
  {
    if (fit_method == "ridge")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$nb_hidden > 0)
      {
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {

            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "lm" || fit_method == "nnls")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ

      if (fit_obj$nb_hidden > 0)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {

            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "mgaussian")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$nb_hidden > 0)
      {
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ
        fit_glmnet <- fit_obj$fit_obj

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w))
          } else {
            newx <- cbind(newx, g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)
            ) %*% w))
          }

          preds <- predict(fit_glmnet,
                           type = "response",
                           s = fit_obj$s,
                           newx = newx)[, , 1]
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        newx <- reformat_cpp(y, lags)
        fit_glmnet <- fit_obj$fit_obj
        preds <- predict(fit_glmnet,
                         type = "response",
                         s = fit_obj$s,
                         newx = newx)[, , 1]
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "glmboost" || fit_method == "xgboost")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      fit_glmboost <- fit_obj$fit_obj
      nb_series <- fit_obj$nb_series
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      scales <- fit_obj$scales
      direct_link <- fit_obj$direct_link

      if (direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- cbind(newx, matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1))
            } else {
              newx <- cbind(newx, matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1)
            } else {
              newx <- matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "gaussian" ||
        fit_method == "matern32" || fit_method == "matern52")
    {
      y <- fit_obj$y
      xm <- fit_obj$xm
      ym <- fit_obj$ym
      xsd <- fit_obj$xsd
      lags <- fit_obj$lags
      xreg <- fit_obj$xreg
      sigma <- fit_obj$sigma
      l <- fit_obj$l
      mat_coefs <- fit_obj$coef
      class_res <- fit_obj$class_res
      series_names <- fit_obj$series_names
      Kxy <- switch(
        fit_method,
        "gaussian" = gaussian_kxy_cpp,
        "matern32" = matern32_kxy_cpp,
        "matern52" = matern52_kxy_cpp
      )

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)
        K_star <- Kxy(
          x = xreg,
          y = my_scale(x = newx, xm = xm,
                       xsd = xsd),
          sigma = sigma,
          l = l
        )
        preds <- drop(crossprod(K_star, mat_coefs)) + ym
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "polynomial")
    {
      y <- fit_obj$y
      xm <- fit_obj$xm
      ym <- fit_obj$ym
      xsd <- fit_obj$xsd
      lags <- fit_obj$lags
      xreg <- fit_obj$xreg
      sigma <- fit_obj$sigma
      l <- fit_obj$l
      d <- fit_obj$d
      mat_coefs <- fit_obj$mat_coefs
      class_res <- fit_obj$class_res
      series_names <- fit_obj$series_names

      for (i in 1:h)
      {
        # enlever l'intercept
        newx <- reformat_cpp(y, lags)
        K_star <-
          poly_kxy_cpp(
            x = xreg,
            y = my_scale(x = newx, xm = xm,
                         xsd = xsd),
            sigma = sigma,
            d = d,
            l = l
          )
        preds <- drop(crossprod(K_star, mat_coefs)) + ym
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "VAR")
    {
      y <- fit_obj$y
      p <- length(fit_obj$series_names)
      if (type_pi == "none")
      {
        preds <- predict(fit_obj$fit_obj, n.ahead = h)
        res <- sapply(1:p,
                      function(i)
                        preds$fcst[[i]][, 1])
        colnames(res) <- fit_obj$series_names
      } else {
        type_pi <- "VAR"
        res <- predict(fit_obj$fit_obj, n.ahead = h,
                       ci = level / 100)
      }
    }

    if (fit_method == "lassoVAR")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      xsd <- fit_obj$scales
      p <- length(fit_obj$series_names)

      for (i in 1:h)
      {
        newinfo <- my_scale(y[1:(lags + 1), ],
                            xm = xm, xsd = xsd)
        newx <- reformat_cpp(newinfo, lags)
        # consider using parSapply when p is high
        preds <- sapply(1:p, function(i) newx%*%fit_obj$fit_lasso[[i]]$beta) + ym
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "scn")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)
        preds <- predict_SCN(fit_obj, newx = newx)
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "pls")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      ym <- fit_obj$ym

      if (fit_obj$direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

             if(fit_obj$hidden_layer_bias == FALSE)
             {
                 newx <- cbind(newx, matrix(g(my_scale(
                   newx, xm = nn_xm,
                   xsd = nn_xsd
                 ) %*% w)[, hidden_layer_index], nrow = 1))
            } else {
                newx <- cbind(newx, matrix(g(my_scale(
                  cbind(1, newx), xm = c(1, nn_xm),
                  xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
            }

          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
          }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1)
          } else {
            newx <- matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
          }

          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "pcr")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      ncomp <- fit_obj$ncomp
      p <- length(fit_obj$series_names)

        for (i in 1:h)
        {
          newx <- my_scale(reformat_cpp(y, lags),
                           xm = fit_obj$xm, xsd = fit_obj$scales)
          preds <- sapply(1:p, function(i) predict(fit_obj$fit_obj[[i]],
                        newdata = newx)[, , ncomp]) + ym
          y <- rbind_vecmat_cpp(preds, y)
        }
    }

    if (fit_method == "gen")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$nb_hidden > 0)
      {
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index],
            nrow = 1))

          } else { # no bias term

            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          cat("\n")
          cat("newx", newx)
          cat("\n")
          print(dim(newx))

          preds <- try(lapply(1:length(fit_obj$series_names),
                              function (i) fit_obj$predict_func(fit_obj$fit_obj[[i]],
                                                                newx = matrix(newx, nrow=1)) + fit_obj$ym[i]),
                                                      silent = TRUE)

          print("here1")

          if (inherits(preds, "try-error"))
          {
            preds <- lapply(1:length(fit_obj$series_names),
                            function (i) fit_obj$predict_func(fit_obj$fit_obj[[i]],
                                                              newdata = matrix(newx, nrow=1)) + fit_obj$ym[i])
          }

          print("here2")

          print(preds)

          y <- rbind_vecmat_cpp(preds, y)
        }
      } else { # nb_hidden <= 0

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- try(fit_obj$predict_func(fit_obj, newx = newx), silent = TRUE)
          if (inherits(preds, "try-error"))
          {
            preds <-fit_obj$predict_func(fit_obj, newdata = newx)
          }
          y <- rbind_vecmat_cpp(preds, y)
        }

      }

    }
  }

  # 2 - direct forecasts -------------------------------------------------

  # direct forecasts
  if (type_forecast == "direct")
  {
    if (fit_method == "ridge2")
    {
      # if col_sample < 1
      if(is.numeric(col_index) == FALSE)
      {
        # observed values (minus lagged)
        y <- fit_obj$y
        # observed values (minus selected columns)
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ


        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_ridge_mts(x = newtrainingx, lags = fit_obj$lags,
                                   nb_hidden = fit_obj$nb_hidden, fit_method = "ridge2",
                                   nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                   hidden_layer_bias = fit_obj$hidden_layer_bias,
                                   col_sample = fit_obj$col_sample, row_sample = fit_obj$row_sample,
                                   a = fit_obj$a, lambda_1 = fit_obj$lambda_1, lambda_2 = fit_obj$lambda_2,
                                   seed = fit_obj$seed)
        }
      } else {
        # observed values (minus lagged)
        y <- fit_obj$y
        # observed values (minus selected columns)
        y_reduced <- as.matrix(fit_obj$y[ , col_index])
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y_reduced, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_ridge_mts(x = newtrainingx, lags = fit_obj$lags,
                                   nb_hidden = fit_obj$nb_hidden, fit_method = "ridge2",
                                   nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                   hidden_layer_bias = fit_obj$hidden_layer_bias,
                                   col_sample = fit_obj$col_sample, row_sample = fit_obj$row_sample,
                                   a = fit_obj$a, lambda_1 = fit_obj$lambda_1, lambda_2 = fit_obj$lambda_2,
                                   seed = fit_obj$seed)
          y_reduced <- as.matrix(y[ , col_index])
        }
      }
    }

    if (fit_method == "lassoVAR")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      xsd <- fit_obj$scales
      lambda <- fit_obj$lambda
      p <- length(fit_obj$series_names)

      for (i in 1:h)
      {
        newinfo <- my_scale(y[1:(lags + 1), ],
                            xm = xm, xsd = xsd)
        newx <- reformat_cpp(newinfo, lags)
        preds <- sapply(1:p, function(i) newx%*%fit_obj$fit_lasso[[i]]$beta) + ym
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
        fit_obj <- fit_var_mts(x = newtrainingx, penalization = "l1",
                               lags = lags, lambda = lambda)
      }
    }

    if (fit_method == "scn")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$method == "greedy")
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_SCN(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_scn_mts(x = newtrainingx, method = "greedy",
                                 lags = fit_obj$lags, activ = fit_obj$activ,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 B = ncol(fit_obj$betas_opt),
                                 nu = fit_obj$nu, lam = fit_obj$lam, r = fit_obj$r,
                                 tol = fit_obj$tol, col_sample = fit_obj$col_sample,
                                 verbose = FALSE,
                                 type_optim = fit_obj$type_optim)
        }
      }

      if (fit_obj$method == "direct")
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_SCN(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_scn_mts(x = newtrainingx, method = "direct",
                                 lags = fit_obj$lags, activ = fit_obj$activ,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 B = ncol(fit_obj$betas_opt),
                                 nu = fit_obj$nu, lam = fit_obj$lam, r = fit_obj$r,
                                 tol = fit_obj$tol, col_sample = fit_obj$col_sample,
                                  verbose = FALSE,
                                 type_optim = fit_obj$type_optim)
        }
      }
    }

    if (fit_method == "pls")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      ym <- fit_obj$ym

      if (fit_obj$direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }
          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_pls_mts(x = newtrainingx,
                                 lags = fit_obj$lags,
                                 activ = fit_obj$activ_name,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 nb_hidden = fit_obj$nb_hidden,
                                 nodes_sim = fit_obj$method,
                                 direct_link = TRUE,
                                 B = fit_obj$ncomp)
        }
      } else { # direct link = FALSE
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1)
          } else {
            newx <- matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
          }

          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_pls_mts(x = newtrainingx,
                                 lags = fit_obj$lags,
                                 activ = fit_obj$activ_name,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 nb_hidden = fit_obj$nb_hidden,
                                 nodes_sim = fit_obj$method,
                                 direct_link = FALSE,
                                 B = fit_obj$ncomp)
        }
      }
    }

    if (fit_method == "pcr")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      ncomp <- fit_obj$ncomp
      p <- length(fit_obj$series_names)

      for (i in 1:h)
      {
        newx <- my_scale(reformat_cpp(y, lags),
                         xm = fit_obj$xm, xsd = fit_obj$scales)
        preds <- sapply(1:p, function(i) predict(fit_obj$fit_obj[[i]],
                                                 newdata = newx)[, , ncomp]) + ym
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1, ]
        fit_obj <- fit_pcr_mts(x = newtrainingx,
                               lags = lags, ncomp = ncomp)
      }
    }

    if (fit_method == "glmboost")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      fit_glmboost <- fit_obj$fit_obj
      nb_series <- fit_obj$nb_series
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      scales <- fit_obj$scales
      direct_link <- fit_obj$direct_link

      if (direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- cbind(newx, matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1))
            } else {
              newx <- cbind(newx, matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          #cat("i = ", i, "\n")
          #cat("preds = ", preds, "\n")
          y <- rbind_vecmat_cpp(preds, y)
          #cat("y = ", y, "\n")
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_glmboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      direct_link = TRUE, a = fit_obj$a, seed = fit_obj$seed)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1)
            } else {
              newx <- matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_glmboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      direct_link = FALSE, a = fit_obj$a, seed = fit_obj$seed)
        }
      }
    }

    if (fit_method == "xgboost")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      fit_xgboost <- fit_obj$fit_obj
      nb_series <- fit_obj$nb_series
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      scales <- fit_obj$scales
      direct_link <- fit_obj$direct_link

      if (direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))

          }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x

          fit_obj <- fit_xgboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                     lambda = fit_obj$lambda, alpha = fit_obj$alpha,
                                     hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      direct_link = TRUE, a = fit_obj$a, seed = fit_obj$seed)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1)
          } else {
            newx <- matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
          }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_xgboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                     hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      direct_link = FALSE, a = fit_obj$a, seed = fit_obj$seed)
        }
      }
    }

    if (fit_method == "gen")
    {
      stop("direct forecast not implemented for 'gen'")
    }
  }

  # Forecast method
  # if no confidence interval is required
  if (type_pi == "none")
  {
    if (fit_method == "VAR") {
      res2 <- rbind(as.matrix(y), res)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n,]
      colnames(res) <- fit_obj$series_names
      return(res)
    } else {
      res2 <- rev_matrix_cpp(y)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n,]
      colnames(res) <- fit_obj$series_names
      return(res)
    }
  } else {
    # if a confidence interval is required
    if (fit_method == "VAR") {
      res2 <- rbind(as.matrix(y),
                    sapply(1:p,
                           function(i)
                             res$fcst[[i]][, 1]))

      ans <- lapply(1:p,
                    function(i) {
                      resid_fcast_matrix <- res$fcst[[i]][, c(1, 2, 3)]
                      colnames(resid_fcast_matrix) <-
                        c("Point Forecast",
                          paste0("Lo ", level),
                          paste0("Hi ", level))
                      resid_fcast_matrix
                    })
      names(ans)  <- fit_obj$series_names

      return(list(obs = res2, fcst = ans))
    } else { # fit_method != "VAR"
      res2 <- rev_matrix_cpp(y)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n,]
      colnames(res) <- fit_obj$series_names
      p <- length(fit_obj$series_names)
      ans <- vector("list", p)
      names(ans) <- fit_obj$series_names
      resid_fcast <- switch(
        type_pi,
        "arima" = lapply(1:p,
                         function(x)
                           forecast::forecast(
                             forecast::auto.arima(fit_obj$resid[, x]),
                             h = h,
                             level = level
                           )),
        "ets" = lapply(1:p,
                       function(x)
                         forecast::forecast(
                           forecast::ets(fit_obj$resid[, x]),
                           h = h,
                           level = level
                         )),
        "theta" = lapply(1:p,
                         function(x)
                           forecast::thetaf(fit_obj$resid[, x],
                                            h = h, level = level))
      )

      nb <- 1 + 2 * length(level)

      for (i in 1:p)
      {
        resid_fcast_matrix <- cbind(0, resid_fcast[[i]]$lower,
                                    resid_fcast[[i]]$upper)
        colnames(resid_fcast_matrix) <-
          c("Point Forecast",
            paste0("Lo ", level),
            paste0("Hi ", level))
        ans[[i]] <- resid_fcast_matrix + matrix(rep(res[, i], nb),
                                                ncol = nb, byrow = FALSE)
      }

      colnames(res2) <- fit_obj$series_names
      return(list(obs = res2, fcst = ans))
    }
  }

}

# regularized regression models
fit_ridge_mts <- function(x,
                          lags = 1,
                          nb_hidden = 5,
                          fit_method = c("ridge", "mgaussian"),
                          nodes_sim = c("sobol", "halton", "unif"),
                          activ = c("relu", "sigmoid", "tanh",
                                    "leakyrelu", "elu", "linear"),
                          hidden_layer_bias = FALSE,
                          col_sample = 1,
                          row_sample = 1,
                          a = 0.01,
                          lambda = 0.1,
                          alpha = 0.5,
                          lambda_1 = 0.1,
                          lambda_2 = 0.1,
                          seed = 1)
{
  stopifnot(col_sample >= 0 && col_sample <= 1)
  stopifnot(floor(col_sample * ncol(x)) >= 1)
  stopifnot(row_sample >= 0.5 && row_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))
  fit_method <- match.arg(fit_method)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  blocks_index <- TRUE
  col_index <- TRUE

  series_names <- colnames(x)

  if (row_sample < 1)
  {
    nb_rows_reduced <- max(4, floor(row_sample * nrow(x)))
    # !!! because the ts object is in reverse order
    set.seed(seed)
    x <- rev_matrix_cpp(x)[1:sample(4:nb_rows_reduced, size = 1),]
  } else {
    # !!! because the ts object is in reverse order
    x <- rev_matrix_cpp(x)
  }

  if (col_sample <= 1)
  {
    nb_cols_reduced <- max(1, floor(col_sample * ncol(x)))
    set.seed(seed)
    col_index <- sort(sample(1:ncol(x), size = nb_cols_reduced))
    blocks_index <- lapply(col_index,
                           function(i)
                             seq(
                               from = (i - 1) * lags + 1,
                               to = i * lags,
                               by = 1
                             ))
    names(blocks_index) <- series_names[col_index]
    blocks_index <- unlist(blocks_index)
  }

  y_x <- create_train_inputs_cpp(x, lags)

  if (nb_hidden > 0)
  {
    list_xreg <- create_new_predictors(
      as.matrix(y_x$xreg[, blocks_index]),
      nb_hidden = nb_hidden,
      hidden_layer_bias = hidden_layer_bias,
      method = nodes_sim,
      activ = activ,
      a = a,
      seed = seed
    )
    xreg <- list_xreg$predictors
  } else {
    xreg <-  y_x$xreg[, TRUE]
  }

  # observed values, minus the lags (!)(beware)(!)
  observed_values <- y <- y_x$y

  if (fit_method == "ridge")
  {
    stopifnot(lambda > 0)
    ym <- colMeans(y)
    centered_y <- my_scale(x = y, xm = ym)
    x_scaled <- my_scale(xreg)
    xreg <- x_scaled$res

    xm <- x_scaled$xm
    xsd <- x_scaled$xsd

    inv <- my_ginv(crossprod(xreg) + lambda * diag(ncol(xreg)))
    lscoef <- inv %*% crossprod(xreg, centered_y)
    colnames(lscoef) <- series_names
    rownames(lscoef) <- colnames(xreg)

    lsfit <- xreg %*% lscoef
    fitted_values <-
      rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                    ncol = ncol(lsfit)))
    resid <- rev_matrix_cpp(observed_values) - fitted_values

    if (nb_hidden > 0)
    {
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = nb_hidden,
          method = nodes_sim,
          w = list_xreg$w,
          activ = list_xreg$activ,
          hidden_layer_bias = hidden_layer_bias,
          hidden_layer_index = list_xreg$hidden_layer_index,
          seed = seed,
          nn_xm = list_xreg$xm,
          nn_xsd = list_xreg$xsd,
          ym = ym,
          xm = xm,
          scales = xsd,
          coef = lscoef,
          resid = resid
        )
      )
    } else {
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = 0,
          ym = ym,
          xm = xm,
          scales = xsd,
          coef = lscoef,
          resid = resid
        )
      )
    }
  }

  if (fit_method == "mgaussian")
  {
    stopifnot(lambda > 0)

    if(is_package_available("glmnet") == FALSE)
    {
      install.packages("glmnet",
                     repos = c(CRAN = "https://cloud.r-project.org"))
    }

    fit_glmnet <-
      glmnet::glmnet(
        x = xreg,
        y = y,
        standardize = FALSE,
        intercept = FALSE,
        alpha = alpha,
        family = "mgaussian"
      )

    fitted_values <-
      predict(fit_glmnet, s = lambda, newx = xreg)[, , 1]

    resid <-
      rev_matrix_cpp(observed_values - fitted_values) # looks false (!)

    if (nb_hidden > 0)
    {
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = nb_hidden,
          method = nodes_sim,
          fit_obj = fit_glmnet,
          s = lambda,
          w = list_xreg$w,
          activ = list_xreg$activ,
          hidden_layer_bias = hidden_layer_bias,
          hidden_layer_index = list_xreg$hidden_layer_index,
          seed = seed,
          nn_xm = list_xreg$xm,
          nn_xsd = list_xreg$xsd,
          resid = resid
        )
      )
    } else {
      # nb_hidden <= 0
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = 0,
          fit_obj = fit_glmnet,
          s = lambda,
          seed = seed,
          resid = resid
        )
      )
    }
  }

}


# least squares regression  model to multiple time series
fit_lm_mts <- function(x,
                       lags = 1,
                       nb_hidden = 5,
                       nodes_sim = c("sobol", "halton", "unif"),
                       activ = c("relu", "sigmoid", "tanh",
                                 "leakyrelu", "elu", "linear"),
                       hidden_layer_bias = FALSE,
                       col_sample = 1,
                       a = 0.01,
                       seed = 1)
{
  stopifnot(col_sample > 0 || col_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  list_xreg <- create_new_predictors(
    y_x$xreg,
    nb_hidden = nb_hidden,
    hidden_layer_bias = hidden_layer_bias,
    method = nodes_sim,
    activ = activ,
    a = a,
    seed = seed
  )
  xreg <- list_xreg$predictors

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  k_p <- lags * ncol(x)
  index <- 1:k_p

  fit_lm <- .lm.fit(x = xreg, y = centered_y)

  lscoef <- fit_lm$coefficients
  colnames(lscoef) <- series_names
  rownames(lscoef) <- colnames(xreg)

  lsfit <- xreg %*% lscoef
  fitted_values <-
    rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                  ncol = ncol(lsfit)))
  resid <- rev_matrix_cpp(observed_values) - fitted_values

  if (nb_hidden > 0)
  {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "lm",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        w = list_xreg$w,
        activ = list_xreg$activ,
        hidden_layer_bias = hidden_layer_bias,
        hidden_layer_index = list_xreg$hidden_layer_index,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  } else {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "lm",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  }
}


# Fitting a constrained regression model
# to multiple time series (with > 0 coeffs)
fit_nnls_mts <- function(x,
                         lags = 1,
                         nb_hidden = 5,
                         nodes_sim = c("sobol", "halton", "unif"),
                         activ = c("relu", "sigmoid", "tanh",
                                   "leakyrelu", "elu", "linear"),
                         hidden_layer_bias = FALSE,
                         col_sample = 1,
                         a = 0.01,
                         seed = 1)
{
  stopifnot(col_sample > 0 || col_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  list_xreg <- create_new_predictors(
    y_x$xreg,
    nb_hidden = nb_hidden,
    method = nodes_sim,
    activ = activ,
    hidden_layer_bias = hidden_layer_bias,
    a = a,
    seed = seed
  )
  xreg <- list_xreg$predictors

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  k_p <- lags * ncol(x)
  index <- 1:k_p

  lscoef <- sapply(1:ncol(centered_y),
                   function (i)
                     coefficients(nnls::nnls(xreg, centered_y[, i])))
  colnames(lscoef) <- series_names
  rownames(lscoef) <- colnames(xreg)

  lsfit <- xreg %*% lscoef
  fitted_values <-
    rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                  ncol = ncol(lsfit)))
  resid <- rev_matrix_cpp(observed_values) - fitted_values

  if (nb_hidden > 0)
  {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "nnls",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        w = list_xreg$w,
        activ = list_xreg$activ,
        hidden_layer_bias = hidden_layer_bias,
        hidden_layer_index = list_xreg$hidden_layer_index,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  } else {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "nnls",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  }
}

# Fitting kernel Ridge regression  model to multiple time series
fit_krls_mts <- function(x,
                         lags = 1,
                         lambda_krls = 0.1,
                         l = 0.1,
                         sigma = 2,
                         d = 1,
                         kernel_type = c("gaussian", "polynomial",
                                         "matern32", "matern52"),
                         inv_method = c("chol", "ginv"))
{
  kernel_type <- match.arg(kernel_type)
  inv_method <- match.arg(inv_method)

  series_names <- colnames(x)

  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  xreg <- y_x$xreg
  x_scaled <- my_scale(xreg)
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd
  y <-  observed_values <- y_x$y
  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)

  n_k <- nrow(xreg)
  p <- ncol(y)

  # Fit with kernel
  K <- switch(
    kernel_type,
    "gaussian" = gaussian_kxx_cpp(x = xreg,
                                  sigma = sigma,
                                  l = l),
    "polynomial" = poly_kxx_cpp(
      x = xreg,
      sigma = sigma,
      d = d,
      l = l
    ),
    "matern32" = matern32_kxx_cpp(x = xreg,
                                  sigma = sigma,
                                  l = l),
    "matern52" = matern52_kxx_cpp(x = xreg,
                                  sigma = sigma,
                                  l = l)
  )

  mat_coefs <- switch(
    inv_method,
    "chol" = chol2inv(chol(K + lambda_krls * diag(n_k))) %*%
      centered_y,
    "ginv" = my_ginv(K + lambda_krls * diag(n_k)) %*% centered_y
  )

  lsfit <- drop(crossprod(K, mat_coefs))
  fitted_values <-
    rev_matrix_cpp(lsfit  + matrix(rep(ym, each = nrow(lsfit)),
                                   ncol = ncol(lsfit)))

  resid <- rev_matrix_cpp(observed_values) - fitted_values

  class_res <- kernel_type

  if (class_res == "gaussian" ||
      class_res == "matern32" || class_res == "matern52")
  {
    return(
      list(
        y = y,
        series_names = series_names,
        lags = lags,
        xreg = xreg,
        sigma = sigma,
        l = l,
        lambda_krls = lambda_krls,
        mat_coefs = mat_coefs,
        class_res = class_res,
        ym = ym,
        xm = xm,
        xsd = xsd,
        resid = resid
      )
    )
  }

  if (class_res == "polynomial")
  {
    return(
      list(
        y = y,
        series_names = series_names,
        lags = lags,
        xreg = xreg,
        sigma = sigma,
        d = d,
        l = l,
        lambda_krls = lambda_krls,
        mat_coefs = mat_coefs,
        class_res = class_res,
        ym = ym,
        xm = xm,
        xsd = xsd
      )
    )
  }
}

# Fitting a VAR model to multiple time series
fit_var_mts <- function(x,
                        lags = 1,
                        penalization = c("none", "l1"),
                        lambda = 0.1,
                        # for penalization == "l1" only
                        type_VAR = c("const", "trend",
                                     "both", "none"))
  # for penalization == "none" only
{
  series_names <- colnames(x)
  lag_names <- as.vector(outer(paste0("lag", 1:lags, "_"),
                               series_names, FUN = "paste0"))
  penalization <- match.arg(penalization)

  # unrestricted VAR algo from package 'vars'
  if (penalization == "none")
  {
    type_VAR <- match.arg(type_VAR)
    fit_obj <- vars::VAR(y = x, p = lags, type = type_VAR)

    resids <- resid(fit_obj$fit_obj)

    return(
      list(
        y = x[-(1:lags), ],
        fit_obj = fit_obj,
        resid = resids,
        series_names = series_names,
        class_res = "VAR"
      )
    )
  }

  # Fu (1998) algo for coordinate descent
  if (penalization == "l1")
  {
    nb_series <- ncol(x)
    # !!! because the ts object is in reverse order
    rev_x <- rev_matrix_cpp(x)

    y_x <- create_train_inputs_cpp(rev_x, 1)
    xreg <- y_x$xreg
    x_scaled <- my_scale(xreg)
    xm <- x_scaled$xm
    xsd <- x_scaled$xsd
    observed_values <- y <- y_x$y
    ym <- colMeans(y)
    centered_y <- my_scale(x = y, xm = ym)

    #X <- reformat_cpp(x_scaled$res,
    #                  n_k = lags)
    X <- x_scaled$res

    XX <- crossprod(X)
    XX2 <- 2 * XX

    # Xy2 <- lapply(1:nb_series,
    #               function (i)
    #                 2 * X * centered_y[1, i])
    Xy2 <- lapply(1:nb_series,
                  function (i)
                    2 * X * centered_y[, i])

    # initial beta parameter for the lasso algo
    # beta0 <- lapply(1:nb_series, function(i) {
    #   coeff <-
    #     as.numeric(solve(XX + lambda * diag(ncol(XX))) %*% (t(X) * centered_y[1, i]))
    #
    #   names(coeff) <- lag_names
    #
    #   return(coeff)
    # })
    beta0 <- lapply(1:nb_series, function(i) {
      coeff <-
        as.numeric(solve(XX + lambda * diag(ncol(XX))) %*% (t(X) * centered_y[, i]))

      names(coeff) <- lag_names

      return(coeff)
    })

    # naming the columns
    names(beta0) <- paste0("y_", series_names)

    # apply lasso algo (shooting)
    # envisage to do this in parallel if there are a lot of series
    fit_lasso <-
      lapply(1:nb_series, function (i)
        lasso_shoot_cpp(
          beta = beta0[[i]],
          XX2 = XX2,
          Xy2 = Xy2[[i]],
          lambda = lambda,
          tol = 1e-05,
          max_iter = 10000
        ))
    names(fit_lasso) <- series_names
    # naming the coefficients
    for (i in 1:nb_series)
      names(fit_lasso[[i]]$beta) <- lag_names

    return(
      list(
        y = y,
        x = as.matrix(x),
        lags = lags,
        lambda = lambda,
        series_names = series_names,
        class_res = "lassoVAR",
        ym = ym,
        xm = xm,
        scales = xsd,
        #resid = NULL,
        fit_lasso = fit_lasso
      )
    )
  }
}

# fitting a pls
fit_pls_mts <- function(x,
                        lags = 1,
                        nb_hidden = 5,
                        nodes_sim = c("sobol", "halton", "unif"),
                        activ = c("relu", "sigmoid", "tanh",
                                  "leakyrelu", "elu", "linear"),
                        hidden_layer_bias = FALSE,
                        direct_link = FALSE,
                        a = 0.01,
                        B = 5,
                        seed = 1)
{
  stopifnot(is.wholenumber(nb_hidden))
  stopifnot(nb_hidden > 0)
  p <- ncol(x)

  series_names <- colnames(x)
  nb_series <- ncol(x)
  # !!! because the ts object is in reverse order
  rev_x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(rev_x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  if (direct_link == FALSE)
  {
    list_xreg <- create_new_predictors(
      as.matrix(y_x$xreg),
      nb_hidden = nb_hidden,
      method = nodes_sim,
      activ = activ,
      hidden_layer_bias = hidden_layer_bias,
      a = a,
      seed = seed
    )
    xreg <- as.matrix(list_xreg$predictors[,-(1:(lags * nb_series))])
  } else {
    list_xreg <- create_new_predictors(
      as.matrix(y_x$xreg),
      nb_hidden = nb_hidden,
      method = nodes_sim,
      activ = activ,
      hidden_layer_bias = hidden_layer_bias,
      a = a,
      seed = seed
    )
    xreg <- list_xreg$predictors
  }

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  ym_mat <- tcrossprod(rep(1, nrow(xreg)), ym)

  fit_obj_pls <- fit_pls(x = xreg, y = centered_y,
                         ncomp = B)
  fitted_values <-
    rev_matrix_cpp(ym_mat + predict_pls(fit_obj = fit_obj_pls,
                                        newx = xreg))
  resids <- fitted_values - rev_matrix_cpp(y)

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      series_names = series_names,
      class_res = "pls",
      nb_hidden = nb_hidden,
      method = nodes_sim,
      w = list_xreg$w,
      activ = list_xreg$activ,
      hidden_layer_bias = hidden_layer_bias,
      hidden_layer_index = list_xreg$hidden_layer_index,
      activ_name = activ,
      direct_link = direct_link,
      seed = seed,
      nn_xm = list_xreg$xm,
      nn_xsd = list_xreg$xsd,
      ym = ym,
      xm = xm,
      scales = xsd,
      coefficients = fit_obj_pls$coefficients,
      Xmeans = fit_obj_pls$Xmeans,
      Ymeans = fit_obj_pls$Ymeans,
      ncomp = B,
      fitted_values = fitted_values,
      resid = resids
    )
  )
}

# fitting a pcr
fit_pcr_mts <- function(x, lags = 1,
                        ncomp = 5)
{
  stopifnot(ncomp <= lags * ncol(x))

  p <- ncol(x)

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  rev_x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(rev_x, lags)
  observed_values <- y <- y_x$y

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(y_x$xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  ym_mat <- tcrossprod(rep(1, nrow(xreg)), ym)

  fit_obj_pcr <- fit_pcr(x = xreg, y = centered_y,
                         ncomp = ncomp)

  fitted_values <- rev_matrix_cpp(ym_mat + sapply(1:p,
                                                  function (i)
                                                    predict(fit_obj_pcr[[i]],
                                                            newdata = xreg)[, , ncomp]))

  resids <- fitted_values - rev_matrix_cpp(y)

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      series_names = series_names,
      class_res = "pcr",
      ym = ym,
      xm = xm,
      scales = xsd,
      ncomp = ncomp,
      fit_obj = fit_obj_pcr,
      fitted_values = fitted_values,
      resid = resids
    )
  )
}

# Check if package is available -----
is_package_available <- function(pkg_name)
{
  return(pkg_name %in% rownames(installed.packages()))
}

# Check if x is an integer -----
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


# Multivariate moving block bootstrap (main loop adapted from Efron and Tibshirani (sec. 8.6)) -----
mbb2 <- function(r,
                 n,
                 b,
                 seed = 123,
                 return_indices = FALSE)
{
  n_obs <- dim(r)[1]
  n_series <- dim(r)[2]
  b <- floor(min(max(3L, b), n_obs - 1L))
  if(n >= n_obs)
    stop("forecasting horizon must be < number of observations")
  n <- min(n_obs, n)

  set.seed(seed) # important for base::sample below

  r_bt <- matrix(NA, nrow = n_obs, ncol = dim(r)[2])  # local vector for a bootstrap replication

  #cat("n_obs", n_obs, "\n")
  #cat("b", b, "\n")
  for (i in 1:ceiling(n_obs/b)) {
    #cat("i: ", i, "----- \n")
    endpoint <- sample(b:n_obs, size = 1)
    #cat("endpoint", endpoint, "\n")
    try(r_bt[(i - 1)*b + 1:b, ] <- r[endpoint - (b:1) + 1, ],
        silent = TRUE)
  }

  tmp <- matrix(r_bt[(1:n), ], nrow = n, ncol = n_series)

  if(return_indices == FALSE)
  {
    return(tmp)
  } else {
    return(arrayInd(match(tmp, r), .dim = dim(r))[1:n, 1])
  }
}


# Multivariate circular block bootstrap (adapted from NMOF book -- Matlab code) -----
mbb <- function(r,
                n,
                b,
                seed = 123,
                return_indices = FALSE)
{
  nT <- dim(r)[1]
  k <- dim(r)[2]

  # b <- (nT + 1)*runif(1); print(b); floor(min(max(2L, b), nT - 1L))
  b <- floor(min(max(3L, b), nT - 1L))

  # circular block bootstrap

  set.seed(seed)

  nb <- ceiling(n / b) # number of bootstrap reps
  js <- floor(runif(n = nb) * nT) # starting points - 1

  x <- matrix(NA, nrow = nb * b, ncol = k)
  for (i in 1:nb)
  {
    j <- ((js[i] + 1:b) %% nT) + 1 #positions in original data
    s <- (1:b) + (i - 1) * b
    x[s, ] <- r[j, ]
  }

  if (nb * n > n)
    # correct length if nb*b > n
  {
    tmp <- drop(x[1:n,])
  } else {
    tmp <- drop(x)
  }

  if (return_indices)
  {
    return(arrayInd(match(tmp, r), .dim = dim(r))[1:nrow(tmp), 1])
  } else {
    return(tmp)
  }
}

#  MASS::ginv -----
my_ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
  {
    stop("'X' must be a numeric or complex matrix")
  }

  Xsvd <- La.svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive))
  {
    return(crossprod(Xsvd$vt, (1 / Xsvd$d * t(Xsvd$u))))
  }
  else if (!any(Positive))
  {
    return(array(0, dim(X)[2L:1L]))
  }
  else {
    return(crossprod(Xsvd$vt[, Positive, drop = FALSE], ((1 / Xsvd$d[Positive]) *
                                                           t(Xsvd$u[, Positive, drop = FALSE]))))
  }
}
my_ginv <- compiler::cmpfun(my_ginv)

# scaling matrices -----
my_scale <- function(x, xm = NULL, xsd = NULL)
{
  rep_1_n <- rep.int(1, dim(x)[1])

  # centering and scaling, returning the means and sd's
  if (is.null(xm) && is.null(xsd))
  {
    xm <- colMeans(x)
    xsd <- my_sd(x)
    return(list(
      res = (x - tcrossprod(rep_1_n, xm)) / tcrossprod(rep_1_n, xsd),
      xm = xm,
      xsd = xsd
    ))
  }

  # centering and scaling
  if (is.numeric(xm) && is.numeric(xsd))
  {
    return((x - tcrossprod(rep_1_n, xm)) / tcrossprod(rep_1_n, xsd))
  }

  # centering only
  if (is.numeric(xm) && is.null(xsd))
  {
    return(x - tcrossprod(rep_1_n, xm))
  }

  # scaling only
  if (is.null(xm) && is.numeric(xsd))
  {
    return(x / tcrossprod(rep_1_n, xsd))
  }
}
my_scale <- compiler::cmpfun(my_scale)

# calculate std's of columns -----
my_sd <- function(x)
{
  n <- dim(x)[1]
  return(drop(rep(1 / (n - 1), n) %*% (x - tcrossprod(
    rep.int(1, n), colMeans(x)
  )) ^ 2) ^ 0.5)
}
my_sd <- compiler::cmpfun(my_sd)

# neutralize -----
neutralize <- function(obj, ym, selected_series)
{
  sims_selected_series <- ahead::getsimulations(obj,
                                                selected_series)$series

  stopifnot(identical(start(sims_selected_series),
                      start(ym)))
  stopifnot(identical(frequency(sims_selected_series),
                      frequency(ym)))
  stopifnot(!is.null(dim(sims_selected_series)))
  stopifnot(is.null(dim(ym)))
  stopifnot(identical(length(ym), nrow(sims_selected_series)))
  stopifnot(all.equal(cumsum(diff(
    time(sims_selected_series)
  )),
  cumsum(diff(time(
    ym
  )))))
  if (any(abs(sims_selected_series) > 2))
    warnings(
      "predictive simulations must contain returns or log-returns \n (use ahead::getreturns on input)"
    )
  n_simulations <- dim(sims_selected_series)[2]
  centered_sims_selected_series <-
    sims_selected_series - matrix(rep(rowMeans(sims_selected_series), n_simulations),
                                  ncol = n_simulations,
                                  byrow = FALSE)
  res <- centered_sims_selected_series + matrix(rep(as.numeric(ym),
                                                    n_simulations),
                                                ncol = n_simulations,
                                                byrow = FALSE)
  return(ts(
    res,
    start = start(sims_selected_series),
    frequency = frequency(sims_selected_series)
  ))
}

# Ridge regression prediction -----
predict_myridge <- function(fit_obj, newx)
{
  my_scale(x = newx,
           xm = fit_obj$xm,
           xsd = fit_obj$scales) %*% fit_obj$coef + fit_obj$ym
}


# Quantile split conformal prediction -----
#'
#' @export
#'
quantile_scp <- function(abs_residuals, alpha)
{
  n_cal_points <- length(abs_residuals)
  k <- ceiling((0.5*n_cal_points + 1)*(1 - alpha))
  return(rank(abs_residuals)[k])
}

# Remove_zero_cols -----
remove_zero_cols <- function(x, with_index = FALSE)
{
  if (with_index == FALSE)
  {
    return(x[, colSums(x == 0) != nrow(x)])
  } else {
    index <- colSums(x == 0) != nrow(x)
    return(list(mat = x[, index],
                index = index))
  }
}

# Fit Ridge regression -----
ridge <- function(x, y, lambda=10^seq(-10, 10,
                                          length.out = 100))
{
  # adapted from MASS::lm.ridge
  x <- as.matrix(x)
  y <- as.vector(y)
  nlambda <- length(lambda)

  ym <- mean(y)
  centered_y <- y - ym

  x_scaled <- base::scale(x)
  attrs <- attributes(x_scaled)
  X <- as.matrix(x_scaled[,])

  Xs <- La.svd(X)
  rhs <- crossprod(Xs$u, centered_y)
  d <- Xs$d
  nb_di <- length(d)
  div <- d ^ 2 + rep(lambda, rep(nb_di, nlambda))
  a <- drop(d * rhs) / div
  dim(a) <- c(nb_di, nlambda)
  n <- nrow(X)

  coef <- crossprod(Xs$vt, a)
  colnames(coef) <- lambda
  centered_y_hat <- X %*% coef

  fitted_values <- drop(ym +  centered_y_hat)
  if (length(lambda) > 1)
  {
    colnames(fitted_values) <- lambda
  }
  residuals <- centered_y - centered_y_hat
  GCV <- colSums(residuals^2)/(n - colSums(matrix(d^2/div, nb_di)))^2
  BIC <- n*log(colMeans(residuals^2)) + (ncol(X) + 2)*log(n)

  out <- list(
      coef = drop(coef),
      ym = ym,
      xm = attrs$`scaled:center`,
      xsd = attrs$`scaled:scale`,
      lambda = lambda,
      best_lam = lambda[which.min(GCV)],
      fitted_values = fitted_values,
      residuals = drop(centered_y - centered_y_hat),
      GCV = GCV,
      BIC = BIC,
      x = x,
      y = y
    )

  return(structure(out, class = "ridge"))
}


predict.ridge <- function(object, newx)
{
  if (length(object$lambda) > 1)
  {
    res <- try(drop(base::scale(newx, center=object$xm,
                                scale=object$xsd)%*%object$coef[,which.min(object$GCV)] + object$ym),
               silent = TRUE)
    if (inherits(res, "try-error"))
    {
      res <- try(drop(base::scale(newx, center=object$xm,
                                  scale=object$xsd)%*%object$coef[which.min(object$GCV)] + object$ym),
                 silent = TRUE)
      return(res)
    } else {
      return(res)
    }
  }  else {
    return(drop(base::scale(newx, center=object$xm,
                            scale=object$xsd)%*%object$coef + object$ym))
  }
}

# Scale a univariate time series -----
scale_ahead <- function(x, center = TRUE, scale = TRUE)
{
  tspx <- tsp(x)
  x <- as.ts(scale.default(x, center = center, scale = scale))
  tsp(x) <- tspx
  return(x)
}

# select residuals distribution -----
select_residuals_dist <- function(resids,
                                  uniformize = c("ranks", "ecdf"),
                                  distro = c("normal", "t", "empirical"))
{
  n_obs <- dim(resids)[1]
  n_series <- dim(resids)[2]
  uniformize <- match.arg(uniformize)
  distro <- match.arg(distro)
  fitted_residuals_distr <- NULL

  if (distro %in% c("normal", "t"))
  {
    fitted_residuals_distr <- vector("list", length = n_series)
    names(fitted_residuals_distr) <- colnames(resids)

    for (j in 1:n_series)
    {
      resid <- resids[,j]
      try_get_res <- suppressWarnings(try(fitdistr_ahead(resid,
                                                         densfun = distro),
                                          silent = TRUE))
      if (!inherits(try_get_res, "try-error"))
      {
        fitted_residuals_distr[[j]] <- try_get_res
      } else {
        warning("distribution can't be fitted by MASS::fitdistr")
        return(NULL)
      }
    }
  }

  if (uniformize == "ranks")
  {
    U <- apply(resids, 2, rank)/(n_obs + 1)
  } else { # uniformize == "ecdf"
    ecdfs <- apply(resids, 2, ecdf)
    U <- sapply(1:n_series, function(i)
      ecdfs[[i]](resids[, i]))
  }

  # select copula between 6 options
  RVM_U <- VineCopula::RVineStructureSelect(data = U, familyset = 1:6,
                                            type = 0, selectioncrit = "BIC")

  return(list(params = fitted_residuals_distr,
              distro = distro,
              RVM_U = RVM_U))
}

# Simulation of residuals using copulas -----
simulate_rvine <- function(obj, RVM_U,
                           h = 5, seed = 123,
                           tests = FALSE)
{
  resids <- obj$resids
  n_obs <- dim(resids)[1]
  n_series <- dim(resids)[2]
  series_names <- colnames(resids)
  params_distro <- obj$params_distro
  distro <- params_distro$distro
  params <- params_distro$params

  # simulate copula
  set.seed(seed)
  rvine_simulation <- VineCopula::RVineSim(N = h, RVM = RVM_U)

  if (identical(distro, "normal"))
  {
    foo <- function (i)
    {
      qnorm(p = rvine_simulation[, i],
            mean = params[[i]]$estimate["mean"],
            sd = params[[i]]$estimate["sd"])
    }

    res <- sapply(1:n_series, function(i) foo(i))
    colnames(res) <- series_names

    if (tests)
    {
      cat("shapiro test p-values: ", "\n")
      tests_p <- apply(res, 2, function(x) stats::shapiro.test(x)$p.value)
      names(tests_p) <- series_names
      print(tests_p)
      cat("\n")
    }
    return(res)
  }

  if (identical(distro, "empirical"))
  {
    foo <- function (i)
    {
      stats::quantile(x = resids[, i],
                      probs = rvine_simulation[, i],
                      type = 7L)
    }

    res <- sapply(1:n_series, function(i) foo(i))
    colnames(res) <- series_names

    return(res)
  }

  if (distro %in% c("student", "t")) {
      stop("Not implemented")
  }

}

# Sort data frame
sort_df <- function(df, by_col)
{
  df[with(df, order(by_col)), ]
}

# Split a time series -----
splitts <- function(y, split_prob = 0.5, return_indices = FALSE, ...)
{
    n_y <- base::ifelse(test = is.null(dim(y)),
                      yes = length(y),
                      no = dim(y)[1])

    index_train <- 1:floor(split_prob*n_y)
    if (return_indices)
      return(index_train)

    start_y <- stats::start(y)
    frequency_y <- stats::frequency(y)

    if(is.null(ncol(y))) # univariate case
    {
        training <- ts(y[index_train],
                       start = start_y,
                       frequency = frequency_y)
        start_testing <- tsp(training )[2] + 1 / frequency_y
        return(list(training = training,
                    testing = ts(y[-index_train],
                                 start = start_testing,
                                 frequency = frequency_y)))
    } else { # multivariate case
      training <- ts(y[index_train, ],
                     start = start_y,
                     frequency = frequency_y)
      start_testing <- tsp(training)[2] + 1 / frequency_y
      return(list(training = training,
                  testing = ts(y[-index_train, ],
                               start = start_testing,
                               frequency = frequency_y)))
    }
}
splitts  <- compiler::cmpfun(splitts)


