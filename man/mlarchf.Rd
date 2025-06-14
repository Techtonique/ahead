% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlarch.R
\name{mlarchf}
\alias{mlarchf}
\title{Conformalized Forecasting using Machine Learning models with ARCH effects}
\usage{
mlarchf(
  y,
  h = 10L,
  mean_model = forecast::auto.arima,
  model_residuals = forecast::thetaf,
  fit_func = ahead::ridge,
  predict_func = predict,
  type_pi = c("kde", "surrogate", "bootstrap"),
  type_sim_conformalize = c("block-bootstrap", "surrogate", "kde", "bootstrap",
    "fitdistr"),
  ml_method = NULL,
  level = 95,
  B = 250L
)
}
\arguments{
\item{y}{A numeric vector or time series of class \code{ts}}

\item{h}{Forecasting horizon}

\item{mean_model}{Function to fit the mean model (default: \code{forecast::auto.arima})}

\item{model_residuals}{Function to model the residuals (default: \code{forecast::thetaf})}

\item{fit_func}{Fitting function for the variance model (default: \code{ahead::ridge})}

\item{predict_func}{Prediction function for the variance model (default: \code{predict})}

\item{type_pi}{Type of prediction interval ("kde", "surrogate", or "bootstrap") for volatility modeling}

\item{type_sim_conformalize}{Type of simulation for conformalization of standardized residuals ("block-bootstrap", "surrogate", "kde", "bootstrap", or "fitdistr")}

\item{ml_method}{Machine learning method to use (if NULL, uses default methods)}

\item{level}{Confidence level for prediction intervals}

\item{B}{Number of bootstrap replications or simulations}
}
\value{
A forecast object containing predictions and prediction intervals
}
\description{
Conformalized Forecasting using Machine Learning models with ARCH effects
}
\examples{
# Example usage will be added
}
