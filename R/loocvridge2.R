#' LOOCV for Ridge2 model
#'
#' LOOCV for Random Vector functional link network model with 2 regularization parameters
#'
#' @param y A multivariate time series of class \code{ts} (preferred) or a \code{matrix}
#' @param xreg External regressors. A data.frame (preferred) or a \code{matrix}
#' @param h Forecasting horizon
#' @param level Confidence level for prediction intervals
#' @param lags Number of lags
#' @param nb_hidden Number of nodes in hidden layer
#' @param nodes_sim Type of simulation for nodes in the hidden layer
#' @param activ Activation function
#' @param a Hyperparameter for activation function "leakyrelu", "elu"
#' @param lambda_1 Regularization parameter for original predictors
#' @param lambda_2 Regularization parameter for transformed predictors
#' @param dropout dropout regularization parameter (dropping nodes in hidden layer)
#' @param seed Reproducibility seed for `nodes_sim == unif`
#' @param type_forecast Recursive or direct forecast
#' @param type_pi Type of prediction interval currently "gaussian", "bootstrap",
#' "blockbootstrap", "movingblockbootstrap", "splitconformal" (very experimental right now),
#' "rvinecopula" (with Gaussian margins for now, Student-t coming soon)
#' @param block_length Length of block for circular or moving block bootstrap
#' @param margins Distribution of margins: "gaussian", "empirical", "student" (postponed or
#' never) for \code{type_pi == "rvinecopula"}
#' @param seed Reproducibility seed for random stuff
#' @param B Number of bootstrap replications or number of simulations (yes, 'B' is unfortunate)
#' @param type_aggregation Type of aggregation, ONLY for bootstrapping; either "mean" or "median"
#' @param centers Number of clusters for \code{type_clustering}
#' @param type_clustering "kmeans" (K-Means clustering) or "hclust" (Hierarchical clustering)
#' @param ym Univariate time series (\code{stats::ts}) of yield to maturities with
#' \code{frequency = frequency(y)} and \code{start = tsp(y)[2] + 1 / frequency(y)}.
#' Default is \code{NULL}.
#' @param cl An integer; the number of clusters for parallel execution, for bootstrap
#' @param show_progress A boolean; show progress bar for bootstrapping? Default is TRUE.
#' @param ... Additional parameters to be passed to \code{\link{kmeans}} or \code{\link{hclust}}
#'
#' @return An object of class "mtsforecast"; a list containing the following elements:
#'
#' \item{method}{The name of the forecasting method as a character string}
#' \item{mean}{Point forecasts for the time series}
#' \item{lower}{Lower bound for prediction interval}
#' \item{upper}{Upper bound for prediction interval}
#' \item{sims}{Model simulations for bootstrapping (basic, or block)}
#' \item{x}{The original time series}
#' \item{residuals}{Residuals from the fitted model}
#' \item{coefficients}{Regression coefficients for \code{type_pi == 'gaussian'} for now}
#'
#' @author T. Moudiki
#'
#' @references
#'
#' Moudiki, T., Planchet, F., & Cousin, A. (2018).
#' Multiple time series forecasting using quasi-randomized
#' functional link neural networks. Risks, 6(1), 22. \cr
#'
#' @export
#'
#' @examples
#'
#' require(fpp)
#'
#' print(ahead::loocvridge2f(fpp::insurance))
#' print(ahead::loocvridge2f(fpp::usconsumption))
#'
#' #foo <- function(xx) ahead::loocvridge2f(fpp::insurance, lambda_1=10^xx[1], lambda_2=10^xx[2])
#' #(opt <- stats::nlminb(objective=foo, lower=c(-10,-10), upper=c(10,10), start=c(0, 0)))
#' #print(ahead::loocvridge2f(fpp::insurance, lambda_1=10^opt$par[1], lambda_2=10^opt$par[2]))
#'
loocvridge2f <- function(y,
                        xreg = NULL,
                        h = 5,
                        level = 95,
                        lags = 1,
                        nb_hidden = 5,
                        nodes_sim = c("sobol", "halton", "unif"),
                        activ = c("relu", "sigmoid", "tanh",
                                "leakyrelu", "elu", "linear"),
                        a = 0.01,
                        lambda_1 = 0.1,
                        lambda_2 = 0.1,
                        dropout = 0,
                        type_forecast = c("recursive", "direct"),
                        type_pi = c(
                        "gaussian",
                        "bootstrap",
                        "blockbootstrap",
                        "movingblockbootstrap",
                        "rvinecopula",
                        "splitconformal"
                        ),
                        block_length = NULL,
                        margins = c("gaussian", "empirical", "student"),
                        seed = 1,
                        B = 100L,
                        type_aggregation = c("mean", "median"),
                        centers = NULL,
                        type_clustering = c("kmeans", "hclust"),
                        ym = NULL,
                        cl = 1L,
                        show_progress = TRUE,
                        ...){
                        nodes_sim <- match.arg(nodes_sim)
                        activ <- match.arg(activ)
                        type_forecast <- match.arg(type_forecast)
                        type_pi <- match.arg(type_pi)
                        margins <- match.arg(margins)
                        type_aggregation <- match.arg(type_aggregation)
                        type_clustering <- match.arg(type_clustering)
                        return(ridge2f(y = y,
                                       xreg = xreg,
                                       h = h,
                                       level = level,
                                       lags = lags,
                                       nb_hidden = nb_hidden,
                                       nodes_sim = nodes_sim,
                                       activ = activ,
                                       a = a,
                                       lambda_1 = lambda_1,
                                       lambda_2 = lambda_2,
                                       dropout = dropout,
                                       type_forecast = type_forecast,
                                       type_pi = type_pi,
                                       block_length = block_length,
                                       margins = margins,
                                       seed = seed,
                                       B = B,
                                       type_aggregation = type_aggregation,
                                       centers = centers,
                                       type_clustering = type_clustering,
                                       ym = ym,
                                       cl = cl,
                                       show_progress = show_progress,
                                       ...))
                    }
