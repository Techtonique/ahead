#' Stacked Doubly-Constrained RVFL for Multivariate Forecasting
#'
#' This function performs a two-stage stacked forecast using ridge regression.
#' In the first stage, it generates base forecasts from the training portion of
#' the time series. In the second stage, these forecasts are used as external
#' regressors to produce the final forecast on the testing set.
#'
#' @param y A multivariate time series object.
#' @param h Integer. Forecast horizon for the second-stage (final) forecast.
#'   Defaults to 5.
#' @param split_fraction Numeric between 0 and 1. Fraction of the time series
#'   used for training. The rest is used for testing and generating stacking
#'   features. Defaults to 0.5.
#' @param ... Additional arguments passed to \code{ahead::ridge2f} in both
#'   stages (e.g., lags, lambda_1, lambda_2, nb_hidden, etc.).
#'
#' @return An object returned by \code{ahead::ridge2f} for the stacked forecast,
#'   typically a list including \code{mean} and prediction intervals.
#'
#' @details
#' The function works as follows:
#' 1. Split the time series into training and testing sets according to
#'    \code{split_fraction}.
#' 2. Generate base forecasts from the training set using \code{ahead::ridge2f}
#' 3. Use these base forecasts as external regressors (\code{xreg}) to predict
#'    the testing set with a second \code{ahead::ridge2f} model.
#' 
#' This approach allows stacking of forecasts to potentially improve accuracy
#' by leveraging the predictions of multiple first-stage models.
#' 
#' For multivariate base learners (e.g., \code{tslm} with multiple dependent
#' variables), the function automatically extracts forecasts from each series
#' and combines them into the feature matrix.
#'
#' @examples
#' \dontrun{
#' 
#' # Univariate example with default ridge2f base learner
#' result1 <- stackridge2f(fpp2::insurance, h = 10, 
#'                         split_fraction = 0.5)
#' }                         
#' @export
stackridge2f <- function(y,
                         h = 5L,
                         split_fraction = 0.5,
                         ...)
{
  split_y <- misc::splitts(y)
  h_testing <- nrow(split_y$testing)
  stage1_forecasts <- ahead::ridge2f(split_y$training, 
                                     h = h_testing)$mean

  return(ahead::ridge2f(
    y = split_y$testing,
    h = h,
    xreg = stage1_forecasts,
    ...))
}