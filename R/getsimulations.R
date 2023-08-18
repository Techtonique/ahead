#' Obtain simulations (when relevant) from a selected time series
#'
#' @param obj result from ridge2f (multivariate time series forecast with simulations)
#' @param selected_series name of the time series selected
#' @param transpose return a transposed time series
#'
#' @export
#'
#' @examples
#'
#' require(fpp)
#'
#' obj <- ahead::ridge2f(fpp::insurance, h = 7,
#'                       type_pi = "bootstrap", B = 5)
#' print(getsimulations(obj, selected_series = "TV.advert"))
#' print(getsimulations(obj, selected_series = "Quotes"))
#' print(getsimulations(obj, selected_series = "TV.advert", transpose = TRUE))
#' print(getsimulations(obj, selected_series = "Quotes", transpose = TRUE))
#'
getsimulations <- function(obj, selected_series, transpose = FALSE)
{
  n_sims <- length(obj$sims)
  start_preds <- stats::start(obj$mean)
  frequency_preds <- stats::frequency(obj$mean)
  selected_series_sims <- sapply(1:n_sims,
                                 FUN = function(i) obj$sims[[i]][, selected_series])

  temp <- ts(selected_series_sims,
             start = start_preds,
             frequency = frequency_preds)

  if (transpose)
  {
    temp <- t(temp)
    colnames(temp) <- paste0("date", 1:nrow(obj$mean))
    return(list(series = temp,
                name = selected_series))
  } else {
    return(list(series = temp,
                name = selected_series))
  }
}
