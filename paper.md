---
title: 'ahead: Univariate and multivariate time series forecasting with uncertainty quantification including simulation approaches'
tags:
- Python
- R
- Julia
- Time series
- Forecasting
- Machine Learning
- Uncertainty
date: "March 25, 2025"
output:
  pdf_document: default
  word_document: default
  html_document:
    df_print: paged
authors:
- name: Thierry Moudiki
  orcid: "0000-0002-9422-5459"
  affiliation: 1
affiliations:
- name: "Techtonique LLC"
  index: 1
bibliography: biblio.bib
---

# Summary

`ahead` contains multiple time series forecasting models, and provides uncertainty quantification through parametric and non-parametric methods, including simulation-based approaches. This paper presents two original Machine Learning models implemented in R version of `ahead` package for forecasting univariate and multivariate time series: `dynrmf`, an autoregressive model that can utilize any Machine Learning model for forecasting univariate time series, and `ridge2f` that extends ridge regression with two regularization parameters and a hidden layer for producing nonlinear outputs.

# Statement of need

Forecasting multivariate time series (MTS hereafter, and including the univariate case) is important for business planning and decision support in finance (stress-testing with financial scenarios), insurance (reserving and required capital valuation), and other industries such as *Energy* (load anticipation) and meteorology. One can obtain point forecasts from a statistical/Machine Learning (ML) model, but these point forecasts are generally of limited importance to analysts. What matters more is the model's ability to quantify the uncertainty around its predictions. 

There are multiple MTS forecasting models available in [\textsf{R} package](https://github.com/Techtonique/ahead) `ahead`'s version `0.18.0` ([\textsf{Python}](https://github.com/Techtonique/ahead_python) and [\textsf{Julia}](https://github.com/Techtonique/Ahead.jl) implementations, following \textsf{R}'s API as closely as possible). R version of `ahead` itself is available at  [r-packages.techtonique.net](https://r-packages.techtonique.net/), and distributed across all major operating systems.

`ahead`'s models generally include parametric prediction intervals alongside non-parametric ones,  simulation-based uncertainty quantification techniques. This paper describes **two** of these ML models, **not available in any other statistical software**:

* \mbox{\texttt{dynrmf}}; an autoregressive dynamic model inspired by **N**eural **N**etwork **A**uto**r**egression (NNAR) (from @hyndman2013forecasting), and described with more details in [https://otexts.com/fpp2/nnetar.html#neural-network-autoregression](https://otexts.com/fpp2/nnetar.html#neural-network-autoregression). As NNAR, \mbox{\texttt{dynrmf}} does an automatic choice of the number of autoregressive and seasonal time series lags. \mbox{\texttt{dynrmf}} is however more generic, and __can use any ML model__ (the default model is [Ridge Regression](https://en.wikipedia.org/wiki/Ridge_regression)).

*  \mbox{\texttt{ahead::ridge2f}} (@moudiki2018multiple) implements a __quasi-randomized *neural* networks__ model extending Ridge_regression to 2 regularization parameters, and capable of producing nonlinear outputs thanks to the use of a *hidden layer*. Since its first publication in 2018, \mbox{\texttt{ahead::ridge2f}} has been enhanced for integrating predictive uncertainty quantification through the (independent/block) bootstrap (@efron1986bootstrap) and copulas' simulation(@brechmann2013modeling, @nagler2023vine). As of march 2025, conformal prediction which improves the quality of prediction intervals (see @vovk2005algorithmic) has been included.

# Examples

## Install `ahead` in \textsf{R}

```r
options(repos = c(
                    techtonique = "https://r-packages.techtonique.net",
                    CRAN = "https://cloud.r-project.org"
                ))
utils::install.packages("rmarkdown", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("remotes", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("forecast", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("fpp", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("ggplot2", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("e1071", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("randomForest", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("dfoptim", repos = c(CRAN="https://cloud.r-project.org"))
utils::install.packages("ahead")
```

```r
library(ahead)
library(forecast)
library(ggplot2)
library(randomForest)
library(e1071)
```

## Use `ahead::ridge2f`

### Use `ahead::ridge2f` for univariate time series forecasting

In all the examples presented below, default hyperparameters are used: 5 nodes in the hidden layer (see @moudiki2018multiple for more details) and a ReLU activation function (see @goodfellow2016deep) are used.

The `fdeaths` data set below contains **monthly deaths of females from various diseases in the UK, 1974-1979**. Here's how to obtain 20-steps-ahead forecasts for `fdeaths` with `ahead::ridge2f`; including seasonality terms. The default level for the prediction interval is equal to 95%.

```r
x <- fdeaths # input dataset
xreg <- ahead::createtrendseason(x) # add seasonality and trend
z <- ahead::ridge2f(x, xreg = xreg, h=20L) # forecasting h-steps ahead
```

```r
ggplot2::autoplot(z) # plot forecast
```
    
![Prediction intervals for `ahead::ridge2f` model on mortality data](paper_files/paper_11_0.png)
    

`EuStocksLogReturns` contains the daily log-returns of major European stock indices, with 1860 observations. Only the first 100 dates of the DAX index dataset are used in the example below.


```r
data(EuStockMarkets)
EuStocks <- ts(EuStockMarkets[1:100, ],
               start = start(EuStockMarkets),
               frequency = frequency(EuStockMarkets)) # original data
EuStocksLogReturns <- ahead::getreturns(EuStocks, type = "log") # obtain log-returns
res <- ahead::ridge2f(EuStocksLogReturns[, "DAX"], h = 20L,
type_pi = "movingblockbootstrap",
                        show_progress = FALSE)
```

```r
ggplot2::autoplot(res) # plot forecast
```
    
![Prediction intervals for `ahead::ridge2f` model on stocks log-returns](paper_files/paper_14_0.png)
    

### Use `ahead::dynrmf` for univariate time series forecasting

`fdeaths` is used in this example too. The default model used by ahead::dynrmf is automated ridge regression, which automatically selects the regularization parameter using Leave-One-Out cross-validation (see @bergmeir2018note):

**- Forecasting with `randomForest::randomForest`**

```r
# Plotting forecasts
# With Random Forest regressor, horizon of 20,
# 95% prediction interval
fit_rf <- ahead::dynrmf(fdeaths, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict)
```

```r
ggplot2::autoplot(fit_rf)
```
    
![Prediction intervals for `ahead::dynrmf` model on mortality data](paper_files/paper_17_0.png)
    

Check in-sample residuals:

```r
forecast::checkresiduals(fit_rf)
```

    
    	Ljung-Box test
    
    data:  Residuals from DynRM 1,1[12]
    Q* = 9.8649, df = 12, p-value = 0.6278
    
    Model df: 0.   Total lags used: 12
    
    
![Residuals diagnostics for `ahead::dynrmf` model on mortality data](paper_files/paper_19_1.png)
    

**- Forecasting with `e1071::svm`** (Support Vector Machines)

```r
# With Support Vector Machine regressor, horizon of 20,
# 95% prediction interval
fit_svm <- ahead::dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
fit_params = list(kernel = "linear"), predict_func = predict)
```

```r
ggplot2::autoplot(fit_svm)
```

![Prediction intervals for `ahead::dynrmf` model on mortality data (using Support Vector Machines regression)](paper_files/paper_22_0.png)
    

Check in-sample residuals:

```r
forecast::checkresiduals(fit_svm)
```

    
    	Ljung-Box test
    
    data:  Residuals from DynRM 1,1[12]
    Q* = 27.351, df = 12, p-value = 0.006875
    
    Model df: 0.   Total lags used: 12
    
![Residuals diagnostics for `ahead::dynrmf` model on mortality data (using Support Vector Machines regression)](paper_files/paper_24_1.png)
    
**- Use of an external regressor** (trend)

The `AirPassengers` dataset has been widely tested as a specialized benchmark in forecasting literature because it displays trend, seasonality, and heteroskedasticity.

```r
h <- 20L # forecasting horizon
res6 <- ahead::dynrmf(AirPassengers, 
xreg_fit = 1:length(AirPassengers),
xreg_predict = (length(AirPassengers)+1):(length(AirPassengers)+h),
h=h)
ggplot2::autoplot(res6)
```
    
![Prediction intervals for `ahead::dynrmf` model on `AirPassengers`](paper_files/paper_26_0.png)
    

### `ahead::ridge2f` for MTS forecasting

The `insurance` dataset (@hyndman2013forecasting) contains monthly quotations and television advertising expenditure for a US insurance company from January 2002 to April 2005. Fast calibration of `ahead::ridge2f` relies on generalized leave-one-out cross-validation as it will be shown in the following \textsf{R} example. It's worth mentioning that **only the 2 regularization parameters are calibrated** here. Other model's hyperparameters such as the number of time series lags or the number of nodes in the hidden layer are set to their default values (respectively `1` and `5`).


```r
objective_function <- function(xx)
{
    ahead::loocvridge2f(fpp::insurance,
                        h = 20L,
                        type_pi="blockbootstrap",
                        lambda_1=10^xx[1],
                        lambda_2=10^xx[2],
                        show_progress = FALSE,
                        )$loocv
}
start <- proc.time()[3]
(opt <- dfoptim::nmkb(fn=objective_function,
                      lower=c(-10,-10),
                      upper=c(10,10),
                      par=c(0.1, 0.1)))
print(proc.time()[3]-start)
```

**Forecasting using the _optimal_ regularization parameters**

The `plot` method (an [S3 method](https://cran.r-project.org/doc/manuals/R-exts.html#Registering-S3-methods)) from `ahead` is used to visualize the predictive simulations and prediction intervals for the `Quotes` and `TV.advert` time series.

```r
start <- proc.time()[3]
res <- ahead::ridge2f(fpp::insurance, h = 20L,
                      type_pi="blockbootstrap",
                      B = 100L, # number of predictive simulations
                      lambda_1=10^opt$par[1], # 'optimal' parameters
                      lambda_2=10^opt$par[2]) # 'optimal' parameters
print(proc.time()[3]-start)


par(mfrow=c(2, 2))
plot(res, "Quotes", type = "sims",
main = "predictive simulations")
plot(res, "TV.advert", type = "sims",
main = "predictive simulations")
plot(res, "Quotes", type = "dist",
main = "prediction intervals")
plot(res, "TV.advert", type = "dist",
main = "prediction intervals")
```
    
![Probabilistic (through simulation) Multivariate forecasting using `ahead::ridge2f`](paper_files/paper_30_1.png)
    

# Citations
