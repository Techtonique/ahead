
debug_print <- function(x) {
  cat("\n")
  print(paste0(deparse(substitute(x)), "'s value:"))
  print(x)
  cat("\n")
}

# simulate Gaussian density -----
rgaussiandens <- function(x,
                          n = length(x),
                          p = 1,
                          seed = 32124) {
  z <- try(stats::density(x, bw = "sj",
                          kernel = "gaussian"),
           silent = TRUE)
  if (inherits(z, "try-error")) {
    z <- density(x, kernel = "gaussian")
  }

  width <- z$bw                              # Kernel width
  rkernel <-
    function(n, seed) {
      set.seed(seed)
      stats::rnorm(n, sd = width)
    }  # Kernel sampler
  if (p <= 1)
  {
    set.seed(seed)
    return(sample(x, n, replace = TRUE) + rkernel(n, seed))    # Here's the entire algorithm
  } else {
    return(sapply(1:p,
                  function(i) {
                    set.seed(seed + i - 1)
                    sample(x, n, replace = TRUE) + rkernel(n, seed + i - 1)
                  }))
  }
}

# split time series -----
splitts <-
  function(y,
           split_prob = 0.5,
           return_indices = FALSE,
           ...)
  {
    n_y <- base::ifelse(test = is.null(dim(y)),
                        yes = length(y),
                        no = dim(y)[1])

    index_train <- 1:floor(split_prob * n_y)
    if (return_indices)
      return(index_train)

    start_y <- stats::start(y)
    frequency_y <- stats::frequency(y)

    if (is.null(ncol(y)))
      # univariate case
    {
      training <- ts(y[index_train],
                     start = start_y,
                     frequency = frequency_y)
      start_testing <- tsp(training)[2] + 1 / frequency_y
      return(list(
        training = training,
        testing = ts(y[-index_train],
                     start = start_testing,
                     frequency = frequency_y)
      ))
    } else {
      # multivariate case
      training <- ts(y[index_train, ],
                     start = start_y,
                     frequency = frequency_y)
      start_testing <- tsp(training)[2] + 1 / frequency_y
      return(list(
        training = training,
        testing = ts(y[-index_train, ],
                     start = start_testing,
                     frequency = frequency_y)
      ))
    }
  }
splitts  <- compiler::cmpfun(splitts)

fit_garch <- function(eps)
{
  return(fGarch::garchFit(
    formula =  ~ garch(1, 1),
    data = eps,
    include.mean = FALSE,
    trace = FALSE
  ))
}

garch11f <- function(y,
                     h = 5,
                     B = 10L,
                     seed = 123)
{
  set.seed(seed)
  if (!is.ts(y))
    y <- ts(y)
  start_y <- stats::start(y)
  frequency_y <- stats::frequency(y)
  tspy <- tsp(y)
  start_preds <- tspy[2] + 1 / tspy[3]
  cat("y: \n", y, "\n")
  fit <- rugarch::ugarchfit(
    data = y,
    spec = rugarch::ugarchspec(),
    solver = "nlminb",
    distribution.model = "norm"
  )
  res <- rugarch::ugarchsim(fit, n.sim = h, m.sim = B)
  print(res@simulation$seriesSim)
  return(list(
    fitted = ts(
      as.numeric(fitted(fit)),
      start = start_y,
      frequency = frequency_y
    ),
    sims = as.matrix(res@simulation$seriesSim)
  ))
}

winkler_score <- function(obj, actual, level = 95) {
  alpha <- 1 - level / 100
  lt <- obj$lower
  ut <- obj$upper
  stopifnot(length(actual) == length(lt) &&
              length(actual) == length(ut))
  n_points <- length(actual)
  diff_lt <- lt - actual
  diff_ut <- actual - ut
  diff_bounds <- ut - lt
  score <- sapply(seq_along(n_points), function(i)
    ifelse(
      actual[i] < lt[i],
      (diff_bounds[i]) + (2 / alpha) * (diff_lt[i]),
      ifelse(
        actual[i] > ut[i],
        (diff_bounds[i]) + (2 / alpha) * (diff_ut[i]),
        # else
        diff_bounds[i]
      )
    ))
  return(mean(score))
}

winkler_score2 <- function(obj, actual, level = 95) {
  alpha <- 1 - level / 100
  lt <- obj$lower
  ut <- obj$upper
  n_points <- length(actual)
  stopifnot((n_points == length(lt)) && (n_points == length(ut)))
  diff_lt <- lt - actual
  diff_bounds <- ut - lt
  diff_ut <- actual - ut
  score <-
    diff_bounds + (2 / alpha) * (pmax(diff_lt, 0) + pmax(diff_ut, 0))
  return(mean(score))
}

coverage_score <- function(obj, actual) {
  if (is.null(obj$lower))
  {
    return(mean((obj$intervals[, 1] <= actual)*(actual <= obj$intervals[, 2]))*100)
  }
  return(mean((obj$lower <= actual)*(actual <= obj$upper))*100)
}

winkler_score <- function(obj, actual, level = 95) {
  alpha <- 1 - level / 100
  lt <- try(obj$lower, silent = TRUE)
  ut <- try(obj$upper, silent = TRUE)
  actual <- as.numeric(actual)
  if (is.null(lt) || is.null(ut))
  {
    lt <- as.numeric(obj$intervals[, 1])
    ut <- as.numeric(obj$intervals[, 2])
  }
  n_points <- length(actual)
  stopifnot((n_points == length(lt)) && (n_points == length(ut)))
  diff_lt <- lt - actual
  diff_bounds <- ut - lt
  diff_ut <- actual - ut
  score <- diff_bounds
  score <- score + (2 / alpha) * (pmax(diff_lt, 0) + pmax(diff_ut, 0))
  return(mean(score))
}

#' Get error metrics
#'
#' This function calculates various error metrics for a given prediction object.
#'
#' @param obj A prediction object containing mean predictions and intervals.
#' @param actual A numeric vector of actual values.
#' @param level The confidence level for the prediction intervals.
#' @return A numeric vector containing the error metrics.
#' @export
#' @examples
#' 
#' obj <- list(mean = 10, lower = 8, upper = 12)
#' actual <- 11
#' geterror(obj, actual)
#' 
geterror <- function(obj, actual, level = 95)
{
  actual <- as.numeric(actual)
  mean_prediction <- as.numeric(obj$mean)
  me <- mean(mean_prediction - actual)
  rmse <- sqrt(mean((mean_prediction - actual)**2))
  mae <- mean(abs(mean_prediction - actual))
  mpe <- mean(mean_prediction/actual-1)
  mape <- mean(abs(mean_prediction/actual-1))
  coverage <- as.numeric(coverage_score(obj, actual))
  winkler <- winkler_score(obj, actual, level = level)
  res <- c(me, rmse, mae, mpe, 
           mape, coverage, winkler)
  names(res) <- c("me", "rmse", "mae", "mpe", 
                  "mape", "coverage", "winkler")
  return(res)
}