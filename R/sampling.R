
#' Simulate using surrogate data
#'
#' @export
#'
rsurrogate <- function(x,
                       n = length(x),
                       p = 1,
                       seed = 123) {
  if (n > length(x))
  {
    stop("For surrogates, must have number of predictions < number of training observations")
  }
  if (p <= 1)
  {
    set.seed(seed)
    res <- tseries::surrogate(x, ns = p,
                              fft = TRUE)
    return(res[seq_len(n), ])
  } else {
    res <- sapply(1:p,
                  function(i) {
                    set.seed(seed + i - 1)
                    tseries::surrogate(x, ns = p,
                                       fft = TRUE)
                  })
    return(res[seq_len(n), ])
  }
}

# fitdistr ------

#' Simulate from parametric distribution
#'
#' @export
#'
rfitdistr <- function(x, n=length(x), p=1)
{
  mean_x <- mean(x)
  sd_x <- sd(x)
  scaled_x <- (x - mean_x)/sd_x
  simulate_function <- misc::fit_param_dist(as.numeric(scaled_x), 
                                            verbose = FALSE)
  res <- simulate_function(n*p)
  if (p <= 1)
  {
    return(res)
  } else {
    return(matrix(res, nrow = n, ncol = p))
  }
}

#' Direct sampling
#' 
#' @param data A numeric vector or matrix.
#' @param n The number of samples to draw.
#' @param method The method to use for sampling.
#' @param kde The kernel density estimate to use for sampling.
#' @param seed The seed to use for sampling.
#' @param ... Additional arguments to pass to the density function.
#'
#' @export
#'
direct_sampling <- function(data = NULL, n = 100L,
                            method = c("kde",
                                       "surrogate",
                                       "bootstrap"),
                            kde = NULL,
                            seed = NULL,
                            ...) {
  method <- match.arg(method)
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  if (identical(method, "kde"))
  {
    if (is.null(kde)) {
      stopifnot(!is.null(data))
      kde <- density(data, bw = "SJ", ...)
    } else if (is.null(data))
    {
      stopifnot(!is.null(kde))
    }
    prob <- kde$y / sum(kde$y)
    return(sample(kde$x, size = n, replace = TRUE, prob = prob))
  }

  if (identical(method, "surrogate"))
  {
    return(sample(tseries::surrogate(data, ns = 1, ...),
                  size = n,
                  replace = TRUE))
  }

  if (identical(method, "bootstrap"))
  {
    return(sample(tseries::tsbootstrap(data, nb = 1, type = "block", b = 1, ...),
                  size = n,
                  replace = TRUE))
  }
}


# Simulate multivariate data -----

#' Simulate multivariate data
#'
#' @param data A numeric vector or matrix.
#' @param method The method to use for sampling.
#' @param n The number of samples to draw.
#' @param block_size The size of the blocks to use for the block bootstrap.
#' @param ... Additional arguments to pass to the density function.
#'
#' @export
#'
rmultivariate <- function(data, method = c("bootstrap", "block-bootstrap"), 
n = 100L, block_size = 5) {
  method <- match.arg(method)
  
  # Ensure data is a matrix
  if (!is.matrix(data)) data <- as.matrix(data)
  n_rows <- nrow(data)
  
  if (method == "bootstrap") {
    # Simple resampling with replacement
    return(data[sample(seq_len(n_rows), size = n, replace = TRUE), ])    
  } 
  
  if (method == "block-bootstrap") {
    # Moving block bootstrap (for time series)
    blocks <- split(1:n_rows, ceiling(seq_along(1:n_rows) / block_size))
    sampled_blocks <- sample(blocks, size = ceiling(n / block_size), replace = TRUE)
    sampled_indices <- unlist(sampled_blocks)[1:n]  # Trim excess if necessary
    return(data[sampled_indices, ])
  }
}
