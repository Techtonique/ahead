# In alphabetical order

# create new predictors
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

# delete columns using a string pattern
delete_columns <- function(x, pattern)
{
  x[,!grepl(pattern = pattern, x = colnames(x))]
}

# dropout regularization
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

# Multivariate circular block bootstrap (adapted from NMOF book -- Matlab code)
mbb <- function(r,
                n,
                b,
                seed = 123,
                return_indices = FALSE)
{
  nT <- dim(r)[1]
  k <- dim(r)[2]

  freq_r <- stats::frequency(r)
  start_r <- stats::start(r)

  stopifnot(b > 1)

  # circular block bootstrap

  set.seed(seed)

  nb <- ceiling(n / b) # number of bootstrap reps
  js <- floor(runif(n = nb) * nT) # starting points - 1

  if (return_indices)
  {
    indices <- vector("list", nb)
    for (i in 1:nb)
    {
      j <- ((js[i] + 1:b) %% nT) + 1 #positions in original data
      indices[[i]] <- j
    }
    return(unlist(indices))
  } else {
    x <- matrix(NA, nrow = nb * b, ncol = k)
    for (i in 1:nb)
    {
      j <- ((js[i] + 1:b) %% nT) + 1 #positions in original data
      s <- (1:b) + (i - 1) * b
      x[s,] <- r[j,]
    }
  }

  if (nb * n > n)
    # correct length if nb*b > n
  {
    return(ts(drop(x[1:n, ]), start = start_r, frequency = freq_r))
  }

  return(ts(drop(x), start = start_r, frequency = freq_r))
}

# adapted from R package forecast (based on Bergmeir et al., 2016, IJF paper)
# mbb2 <- function (x, window_size, seed = 123, return_indices = FALSE)
# {
#   stopifnot(!is.null(ncol(x)))
#   n_series <- dim(x)[2]
#   n_obs_series <- dim(x)[1]
#   set.seed(seed)
#   n_buckets_indices <- floor(n_obs_series / window_size) + 2
#   bx <- matrix(0, nrow = n_buckets_indices * window_size,
#                ncol = n_series)
#   if (return_indices)
#   {
#     indices <- vector("list", n_buckets_indices)
#
#     for (i in 1:n_buckets_indices) {
#       c_ <- sample(1:(n_obs_series - window_size + 1), 1)
#       j <- c_:(c_ + window_size - 1) #positions in original data
#       indices[[i]] <- j # this doesn't work (one more step's needed)
#     }
#
#     return(unlist(indices))
#
#   } else {
#
#     for (i in 1:n_buckets_indices) {
#       c_ <- sample(1:(n_obs_series - window_size + 1), 1)
#       bx[((i - 1) * window_size + 1):(i * window_size), ] <- x[c_:(c_ +
#                                                                     window_size - 1), ]
#     }
#   }
#
#   start_from <- sample(0:(window_size - 1), 1) + 1
#   return(bx[start_from:(start_from + n_obs_series - 1), ])
# }

#  MASS::ginv
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


# scaling matrices
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


# calculate std's of columns
my_sd <- function(x)
{
  n <- dim(x)[1]
  return(drop(rep(1 / (n - 1), n) %*% (x - tcrossprod(
    rep.int(1, n), colMeans(x)
  )) ^ 2) ^ 0.5)
}
my_sd <- compiler::cmpfun(my_sd)


# Ridge regression prediction
predict_myridge <- function(fit_obj, newx)
{
  my_scale(x = newx,
           xm = fit_obj$xm,
           xsd = fit_obj$scales) %*% fit_obj$coef + fit_obj$ym
}


# Remove_zero_cols
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


# Scale a univariate time series
scale_ahead <- function(x, center = TRUE, scale = TRUE) {
  tspx <- tsp(x)
  x <- as.ts(scale.default(x, center = center, scale = scale))
  tsp(x) <- tspx
  return(x)
}
