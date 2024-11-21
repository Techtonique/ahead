direct_sampling <- function(data = NULL, n = 1000,
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