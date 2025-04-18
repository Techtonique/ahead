% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridge2.R
\name{ridge2f}
\alias{ridge2f}
\title{Ridge2 model}
\usage{
ridge2f(
  y,
  h = 5,
  level = 95,
  xreg = NULL,
  lags = 1,
  nb_hidden = 5,
  nodes_sim = c("sobol", "halton", "unif"),
  activ = c("relu", "sigmoid", "tanh", "leakyrelu", "elu", "linear"),
  a = 0.01,
  lambda_1 = 0.1,
  lambda_2 = 0.1,
  dropout = 0,
  type_forecast = c("recursive", "direct"),
  type_pi = c("gaussian", "bootstrap", "blockbootstrap", "movingblockbootstrap",
    "rvinecopula", "conformal-split", "conformal-bootstrap", "conformal-block-bootstrap"),
  block_length = NULL,
  margins = c("gaussian", "empirical"),
  seed = 1,
  B = 100L,
  type_aggregation = c("mean", "median"),
  centers = NULL,
  type_clustering = c("kmeans", "hclust"),
  ym = NULL,
  cl = 1L,
  show_progress = TRUE,
  ...
)
}
\arguments{
\item{y}{A univariate of multivariate time series of class \code{ts} (preferred) or a \code{matrix}}

\item{h}{Forecasting horizon}

\item{level}{Confidence level for prediction intervals}

\item{xreg}{External regressors. A data.frame (preferred) or a \code{matrix}}

\item{lags}{Number of lags}

\item{nb_hidden}{Number of nodes in hidden layer}

\item{nodes_sim}{Type of simulation for nodes in the hidden layer}

\item{activ}{Activation function}

\item{a}{Hyperparameter for activation function "leakyrelu", "elu"}

\item{lambda_1}{Regularization parameter for original predictors}

\item{lambda_2}{Regularization parameter for transformed predictors}

\item{dropout}{dropout regularization parameter (dropping nodes in hidden layer)}

\item{type_forecast}{Recursive or direct forecast}

\item{type_pi}{Type of prediction interval currently "gaussian", "bootstrap",
"blockbootstrap", "movingblockbootstrap", "rvinecopula" (with Gaussian and empirical margins for now)
"conformal-split", "conformal-bootstrap", "conformal-block-bootstrap"}

\item{block_length}{Length of block for circular or moving block bootstrap}

\item{margins}{Distribution of margins: "gaussian", "empirical" for \code{type_pi == "rvinecopula"}}

\item{seed}{Reproducibility seed for random stuff}

\item{B}{Number of bootstrap replications or number of simulations (yes, 'B' is unfortunate)}

\item{type_aggregation}{Type of aggregation, ONLY for bootstrapping; either "mean" or "median"}

\item{centers}{Number of clusters for \code{type_clustering}}

\item{type_clustering}{"kmeans" (K-Means clustering) or "hclust" (Hierarchical clustering)}

\item{ym}{Univariate time series (\code{stats::ts}) of yield to maturities with
\code{frequency = frequency(y)} and \code{start = tsp(y)[2] + 1 / frequency(y)}.
Default is \code{NULL}.}

\item{cl}{An integer; the number of clusters for parallel execution, for bootstrap}

\item{show_progress}{A boolean; show progress bar for bootstrapping? Default is TRUE.}

\item{...}{Additional parameters to be passed to \code{\link{kmeans}} or \code{\link{hclust}}}
}
\value{
An object of class "mtsforecast"; a list containing the following elements:

\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts for the time series}
\item{lower}{Lower bound for prediction interval}
\item{upper}{Upper bound for prediction interval}
\item{sims}{Model simulations for bootstrapping (basic, or block)}
\item{x}{The original time series}
\item{residuals}{Residuals from the fitted model}
\item{coefficients}{Regression coefficients for \code{type_pi == 'gaussian'} for now}
}
\description{
Random Vector functional link network model with 2 regularization parameters
}
\examples{

require(fpp)

print(ahead::ridge2f(fpp::insurance)$mean)
print(ahead::ridge2f(fpp::usconsumption)$lower)

res <- ahead::ridge2f(fpp::insurance, h=10, lags=2)
par(mfrow=c(1, 2))
plot(res, "Quotes")
plot(res, "TV.advert")

# include a trend (just for the example)
xreg <- as.numeric(time(fpp::insurance))
res2 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
h=10, lags=2)
par(mfrow=c(1, 2))
plot(res2, "Quotes")
plot(res2, "TV.advert")

# block bootstrap
xreg <- as.numeric(time(fpp::insurance))
res3 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
                      h=10, lags=1L, type_pi = "bootstrap", B=10)
res5 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
                      h=10, lags=1L, type_pi = "blockbootstrap", B=10,
                      block_length = 4)

print(res3$sims[[2]])
print(res5$sims[[2]])

par(mfrow=c(2, 2))
plot(res3, "Quotes")
plot(res3, "TV.advert")
plot(res5, "Quotes")
plot(res5, "TV.advert")


res4 <- ahead::ridge2f(fpp::usconsumption, h=20, lags=2L,
lambda_2=1)
par(mfrow=c(1, 2))
plot(res4, "income")
plot(res4, "consumption")


# moving block bootstrap
xreg <- as.numeric(time(fpp::insurance))
res6 <- ahead::ridge2f(fpp::insurance, xreg=xreg,
                      h=10, lags=1L,
                      type_pi = "movingblockbootstrap", B=10,
                      block_length = 4)

print(res6$sims[[2]])

par(mfrow=c(1, 2))
plot(res6, "Quotes")
plot(res6, "TV.advert")


}
\references{
Moudiki, T., Planchet, F., & Cousin, A. (2018).
Multiple time series forecasting using quasi-randomized
functional link neural networks. Risks, 6(1), 22. \cr
}
\author{
T. Moudiki
}
