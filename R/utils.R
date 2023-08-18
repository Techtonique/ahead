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
  b <- floor(min(max(2L, b), nT - 1L))

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

# Split a time series -----
splitts <- function(y, p = 0.8, return_indices = TRUE, ...)
{
    n_y <- base::ifelse(test = is.null(dim(y)),
                      yes = length(y),
                      no = dim(y)[1])

    index_train <- 1:floor(p*n_y)
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

