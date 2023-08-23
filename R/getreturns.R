#' Calculate returns or log-returns for multivariate time series
#'
#' @param x Multivariate time series
#' @param type Type of return: basic return ("basic") or log-return ("log")
#'
#' @return The returns
#' @export
#'
#' @examples
#'
#' returns <- getreturns(EuStockMarkets)
#' log_returns <- getreturns(EuStockMarkets,
#'                           type = "log")
#'
#' par(mfrow=c(1, 3))
#' matplot(EuStockMarkets, type = 'l', main = "Closing Prices of \n European stocks (1991-1998)",
#' xlab = "time")
#' matplot(returns, type = 'l', main = "Returns", xlab = "time")
#' matplot(log_returns, type = 'l', main = "Log-returns", xlab = "time")
#'
getreturns <- function(x, type = c("basic", "log"))
{
  # ?rlang::abort #?
  stopifnot(is.ts(x))
  names_x <- colnames(x)
  type <- match.arg(type)

  if (identical(type, "basic"))
  {
    res <- diff(x)/lag(x)
    colnames(res) <- names_x
    return (res)
  }
  res <- diff(log(x))
  colnames(res) <- names_x
  return (res)
}
