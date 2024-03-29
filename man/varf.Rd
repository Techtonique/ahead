% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VAR.R
\name{varf}
\alias{varf}
\title{Vector Autoregressive model (adapted from vars::VAR)}
\usage{
varf(
  y,
  h = 5,
  level = 95,
  lags = 1,
  type_VAR = c("const", "trend", "both", "none"),
  ...
)
}
\arguments{
\item{y}{A multivariate time series of class \code{ts}}

\item{h}{Forecasting horizon}

\item{level}{Confidence level for prediction intervals}

\item{lags}{Number of lags}

\item{type_VAR}{Type of deterministic regressors to include.}

\item{...}{Additional parameters to be passed to vars::VAR.}
}
\value{
An object of class "mtsforecast"; a list containing the following elements:

\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts for the time series}
\item{lower}{Lower bound for prediction interval}
\item{upper}{Upper bound for prediction interval}
\item{x}{The original time series}
\item{residuals}{Residuals from the fitted model}
}
\description{
Vector Autoregressive model adapted from vars::VAR (only for benchmarking)
}
\examples{

require(fpp)

print(varf(fpp::insurance, lags=2, h=10))

res <- varf(fpp::usconsumption, h=20, lags=2)

par(mfrow=c(1, 2))
plot(res, "consumption")
plot(res, "income")

}
\references{
Bernhard Pfaff (2008). VAR, SVAR and SVEC Models: Implementation
Within R Package vars. Journal of Statistical Software 27(4). URL
http://www.jstatsoft.org/v27/i04/. \cr

Pfaff, B. (2008) Analysis of Integrated and Cointegrated Time Series
with R. Second Edition. Springer, New York. ISBN 0-387-27960-1
}
\author{
T. Moudiki
}
