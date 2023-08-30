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
    return(arrayInd(match(tmp, r), .dim = dim(r))[1:dim(tmp)[1], 1])
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
    ecdfs <- apply(obj$residuals, 2, ecdf)
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


# Stratify stuff -----
# from https://gist.github.com/mrdwab/6424112
stratified <- function(df, group, size, select = NULL,
                       replace = FALSE, bothSets = FALSE) {
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      # message(
      #   "Some groups\n---",
      #   paste(names(df.table[df.table < size]), collapse = ", "),
      #   "---\ncontain fewer observations",
      #   " than desired number of samples.\n",
      #   "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)

  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}


