# In alphabetical order

# compute prediction intervals for univariate time series -----

# compute_uni_pi <- function(out, h = 5,
#                            type_pi = c("E", "A", "T", "gaussian"))
# {
#   tspx <- tsp(out$x)
#
#   if (type_pi == "E")
#   {
#     resid_fcast <- forecast::forecast(forecast::ets(out$residuals),
#                                       h = h, level=out$level)
#   }
#
#   if (type_pi == "A")
#   {
#     resid_fcast <- forecast::forecast(forecast::auto.arima(out$residuals),
#                                       h = h, level=out$level)
#   }
#
#   if (type_pi == "T")
#   {
#     resid_fcast <- forecast::thetaf(out$residuals, h = h, level=out$level)
#   }
#
#   if (type_pi == "gaussian")
#   {
#     qt_sd <- -qnorm(0.5 - out$level/200)*sd(out$residuals)
#     rep_qt_sd <- ts(rep(qt_sd, h), start = tspx[2] + 1 / tspx[3],
#                     frequency = tspx[3])
#
#     resid_fcast <- list()
#     resid_fcast$mean <- 0
#     resid_fcast$lower <- -rep_qt_sd
#     resid_fcast$upper <- rep_qt_sd
#   }
#
#   out$mean <-  drop(out$mean  + resid_fcast$mean)
#   out$lower <- drop(out$mean + resid_fcast$lower)
#   out$upper <- drop(out$mean + resid_fcast$upper)
#
#   return(structure(out, class = "forecast"))
# }

# create new predictors -----
create_new_predictors <- function(x,
                                  nb_hidden = 5,
                                  method = c("sobol", "halton", "unif"),
                                  activ = c("relu", "sigmoid", "tanh",
                                            "leakyrelu", "elu", "linear"),
                                  a = 0.01,
                                  seed = 123)
{
  n <- nrow(x)
  stopifnot(nb_hidden > 0)

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

# prehistoric stuff -----
debug_print <- function(x) {
  cat("\n")
  print(paste0(deparse(substitute(x)), "'s value:"))
  print(x)
  cat("\n")
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
      stop("uncomment")
      # if (!is.null(start))
      #   stop(gettextf("supplying pars for the %s distribution is not supported",
      #                 "log-Normal"), domain = NA)
      # if (any(x <= 0))
      #   stop("need positive values to fit a log-Normal")
      # lx <- log(x)
      # sd0 <- sqrt((n - 1)/n) * sd(lx)
      # mx <- mean(lx)
      # estimate <- c(mx, sd0)
      # sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      # names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      # vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2,
      #              dimnames = list(names(sds), names(sds)))
      # names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      # return(structure(list(estimate = estimate, sd = sds,
      #                       vcov = vc, n = n, loglik = sum(stats::dlnorm(x, mx,
      #                                                             sd0, log = TRUE))), class = "fitdistr"))
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
    # if (distname == "poisson") {
    #   if (!is.null(start))
    #     stop(gettextf("supplying pars for the %s distribution is not supported",
    #                   "Poisson"), domain = NA)
    #   estimate <- mean(x)
    #   sds <- sqrt(estimate/n)
    #   names(estimate) <- names(sds) <- "lambda"
    #   vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("lambda",
    #                                                           "lambda"))
    #   return(structure(list(estimate = estimate, sd = sds,
    #                         vcov = vc, n = n, loglik = sum(stats::dpois(x, estimate,
    #                                                              log = TRUE))), class = "fitdistr"))
    # }
    # if (distname == "exponential") {
    #   if (any(x < 0))
    #     stop("Exponential values must be >= 0")
    #   if (!is.null(start))
    #     stop(gettextf("supplying pars for the %s distribution is not supported",
    #                   "exponential"), domain = NA)
    #   estimate <- 1/mean(x)
    #   sds <- estimate/sqrt(n)
    #   vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("rate",
    #                                                           "rate"))
    #   names(estimate) <- names(sds) <- "rate"
    #   return(structure(list(estimate = estimate, sd = sds,
    #                         vcov = vc, n = n, loglik = sum(stats::dexp(x, estimate,
    #                                                             log = TRUE))), class = "fitdistr"))
    # }
    # if (distname == "geometric") {
    #   if (!is.null(start))
    #     stop(gettextf("supplying pars for the %s distribution is not supported",
    #                   "geometric"), domain = NA)
    #   estimate <- 1/(1 + mean(x))
    #   sds <- estimate * sqrt((1 - estimate)/n)
    #   vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("prob",
    #                                                           "prob"))
    #   names(estimate) <- names(sds) <- "prob"
    #   return(structure(list(estimate = estimate, sd = sds,
    #                         vcov = vc, n = n, loglik = sum(stats::dgeom(x, estimate,
    #                                                              log = TRUE))), class = "fitdistr"))
    # }
    # if (distname == "weibull" && is.null(start)) {
    #   if (any(x <= 0))
    #     stop("Weibull values must be > 0")
    #   lx <- log(x)
    #   m <- mean(lx)
    #   v <- stats::var(lx)
    #   shape <- 1.2/sqrt(v)
    #   scale <- exp(m + 0.572/shape)
    #   start <- list(shape = shape, scale = scale)
    #   start <- start[!is.element(names(start), dots)]
    # }
    # if (distname == "gamma" && is.null(start)) {
    #   if (any(x < 0))
    #     stop("gamma values must be >= 0")
    #   m <- mean(x)
    #   v <- stats::var(x)
    #   start <- list(shape = m^2/v, rate = m/v)
    #   start <- start[!is.element(names(start), dots)]
    #   control <- list(parscale = c(1, start$rate))
    # }
    # if (distname == "negative binomial" && is.null(start)) {
    #   m <- mean(x)
    #   v <- stats::var(x)
    #   size <- if (v > m)
    #     m^2/(v - m)
    #   else 100
    #   start <- list(size = size, mu = m)
    #   start <- start[!is.element(names(start), dots)]
    # }
    # if (is.element(distname, c("cauchy", "logistic")) &&
    #     is.null(start)) {
    #   start <- list(location = median(x), scale = stats::IQR(x)/2)
    #   start <- start[!is.element(names(start), dots)]
    # }
    if (distname == "t" && is.null(start)) {
      start <- list(m = median(x), s = stats::IQR(x)/2, df = 10)
      start <- start[!is.element(names(start), dots)]
    }
  }
  # if (is.null(start) || !is.list(start))
  #   stop("'start' must be a named list")
  # nm <- names(start)
  # f <- formals(densfun)
  # args <- names(f)
  # m <- match(nm, args)
  # if (any(is.na(m)))
  #   stop("'start' specifies names which are not arguments to 'densfun'")
  # formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
  # dens <- function(parm, x, ...) densfun(x, parm, ...)
  # if ((l <- length(nm)) > 1L)
  #   body(dens) <- parse(text = paste("densfun(x,", paste("parm[",
  #                                                        1L:l, "]", collapse = ", "), ", ...)"))
  #
  # Call[[1L]] <- quote(stats::nlminb)
  # Call$densfun <- Call$start <- NULL
  # Call$x <- x
  # Call$start <- start
  #
  # Call$objective <- if ("log" %in% args)
  #   mylogfn
  # else myfn
  #
  # Call$hessian <- TRUE # in MASS::fitdistr
  # # Call$hessian <- NULL
  # if (length(control))
  #   Call$control <- control
  #
  # res <- eval.parent(Call)
  #
  # return(list(estimate = res$par,
  #             convergence = res$convergence,
  #             objective = res$objective))
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

# Check if package is available -----
is_package_available <- function(pkg_name)
{
  return(pkg_name %in% rownames(installed.packages()))
}

# Check if x is an integer -----
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


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

  tmp <- drop(x[1:n,])

  if (return_indices)
  {
    return(arrayInd(match(tmp, r), .dim = dim(r))[1:nrow(tmp), 1])
  } else {
    return(tmp)
  }
}


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
#' @export 
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

#' @export
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

#' @export
removenas <- function(y) {
  # Check if input is a time series object
  if (!is.ts(y)) {
    stop("Input must be a time series object (ts).")
  }

  # Check for NA values
  if (!anyNA(y)) {
    return(y)
  }

  # Extract time and values
  time <- time(y)
  values <- as.numeric(y)

  # Perform linear interpolation to replace NAs
  interpolated_values <- approx(x = time[!is.na(values)], y = values[!is.na(values)], xout = time)$y

  # Create a new time series object
  new_ts <- ts(interpolated_values, start = start(y), frequency = frequency(y))

  return(new_ts)
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
                                  distro = c("normal", "empirical"))
{
  n_obs <- dim(resids)[1]
  n_series <- dim(resids)[2]
  uniformize <- match.arg(uniformize)
  distro <- match.arg(distro)
  fitted_residuals_distr <- NULL

  if (identical(distro, "normal"))
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

  # if (distro %in% c("student", "t")) {
  #     stop("Not implemented")
  # }

}

# Sort data frame
sort_df <- function(df, by_col)
{
  df[with(df, order(by_col)), ]
}

# Split a time series -----
# splitts <- function(y, split_prob = 0.5, return_indices = FALSE, ...)
# {
#     n_y <- base::ifelse(test = is.null(dim(y)),
#                       yes = length(y),
#                       no = dim(y)[1])
#
#     index_train <- 1:floor(split_prob*n_y)
#     if (return_indices)
#       return(index_train)
#
#     start_y <- stats::start(y)
#     frequency_y <- stats::frequency(y)
#
#     if(is.null(ncol(y))) # univariate case
#     {
#         training <- ts(y[index_train],
#                        start = start_y,
#                        frequency = frequency_y)
#         start_testing <- tsp(training )[2] + 1 / frequency_y
#         return(list(training = training,
#                     testing = ts(y[-index_train],
#                                  start = start_testing,
#                                  frequency = frequency_y)))
#     } else { # multivariate case
#       training <- ts(y[index_train, ],
#                      start = start_y,
#                      frequency = frequency_y)
#       start_testing <- tsp(training)[2] + 1 / frequency_y
#       return(list(training = training,
#                   testing = ts(y[-index_train, ],
#                                start = start_testing,
#                                frequency = frequency_y)))
#     }
# }
# splitts  <- compiler::cmpfun(splitts)
#
#
