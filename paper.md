# 1 - Introduction and context

Univariate and Multivariate time series (MTS hereafter) are collections
of sequential data points observed at different timesteps for measurable
indicators. Real-world examples of MTS include – among many others – the
monthly totals of international airline passengers from 1949 to 1960, or
the daily closing prices of major European stock indices from 1991 to
1998.

Forecasting MTS is important for business planning and decision support
in finance, insurance, and other industries such as *Energy* (load
anticipation) and meteorology. One can obtain point forecasts from a
statistical/Machine Learning model, but these point forecasts are
generally of limited importance to analysts. What matters more is the
model’s ability to quantify the uncertainty around its predictions. In
finance or insurance for example, uncertainty quantification through
simulation is crucial for risk management, liabilities’ reserving, and
capital valuation.

There are multiple MTS forecasting methods implemented in [
package](https://github.com/Techtonique/ahead) `ahead`’s version
`0.11.0` (there are also [](https://github.com/Techtonique/ahead_python)
and [](https://github.com/Techtonique/Ahead.jl) implementations, both
following ’s API as closely as possible). `ahead` is available through
the [R-universe](https://techtonique.r-universe.dev/builds), which
allows the package to be continuously integrated and distributed across
all major platforms (Windows, macOS, Linux).

All of `ahead`’s forecasting methods include parametric prediction
intervals alongside non-parametric, simulation-based uncertainty
quantification techniques. This paper describes **two** of these Machine
Learning methods, the *original* ones, **not available in any other
statistical software**:

-   ; an autoregressive dynamic Machine Learning model inspired by NNAR
    (**N**eural **N**etwork **A**uto**r**egression, see Hyndman and
    Athanasopoulos (2013). As NNAR, does an automatic choice of the
    number of autoregressive and seasonal time series lags. However,
    instead of an artificial neural network, as implemented in NNAR,
    **can use any regression (Machine Learning) model**.

-   The forecasting model from Moudiki, Planchet, and Cousin (2018): .
    implements a **quasi-randomized *neural* networks** model extending
    [ridge regression](https://en.wikipedia.org/wiki/Ridge_regression)
    with 2 regularization parameters, and capable of producing highly
    nonlinear outputs through the use of a *hidden layer*. Since its
    first publication in 2018, has been enhanced for integrating
    uncertainty quantification through the (independent, circular,
    moving block) bootstrap (Efron and Tibshirani (1986)) and copulas
    simulations (Brechmann and Schepsmeier (2013), Nagler et al.
    (2023)). Future (ongoing) developments include conformal prediction
    (Vovk, Gammerman, and Shafer (2005)) and Kernel Density Estimation
    (Silverman (2018)).

# 2 - Examples of use

## 2 - 1 How to install `ahead` in 

    options(repos = c(
        techtonique = 'https://techtonique.r-universe.dev',
        CRAN = 'https://cloud.r-project.org'))
    utils::install.packages("rmarkdown", repos = c(CRAN="https://cloud.r-project.org"))
    utils::install.packages("remotes", repos = c(CRAN="https://cloud.r-project.org"))
    utils::install.packages("forecast", repos = c(CRAN="https://cloud.r-project.org"))
    utils::install.packages("ggplot2", repos = c(CRAN="https://cloud.r-project.org"))
    utils::install.packages("e1071", repos = c(CRAN="https://cloud.r-project.org"))
    utils::install.packages("randomForest", repos = c(CRAN="https://cloud.r-project.org"))
    utils::install.packages(c("ahead", "dfoptim"))

**Loading packages**

    suppressPackageStartupMessages(library(ahead))
    suppressPackageStartupMessages(library(forecast))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(randomForest))
    suppressPackageStartupMessages(library(e1071))

## 2 - 2 Use `ahead::ridge2f`

## 2 - 2 - 1 Use `ahead::ridge2f` for univariate time series forecasting

In all these examples, in order to capture the nonlinear patterns of the
inputs, 5 nodes in the hidden layer and a ReLU activation function are
used (default hyperparameters, see Goodfellow, Bengio, and Courville
(2016) and Moudiki, Planchet, and Cousin (2018) for more details).

The `fdeaths` data set below contains **monthly deaths of females from
bronchitis, emphysema and asthma in the UK, 1974-1979**. The data were
collected by the Office of Population Censuses and Surveys. Mortality
data are widely used in insurance, actuarial science and demography.
They are notably useful for the valuation of death-linked liabilities.

Here’s how to obtain 20-steps-ahead forecasts for `fdeaths` with
`ahead::ridge2f`; including seasonality terms. The default level for the
prediction interval is 95%.

    x <- fdeaths # input dataset
    xreg <- ahead::createtrendseason(x) # add seasonality and trend
    z <- ahead::ridge2f(x, xreg = xreg, h=20L) # forecasting h-steps ahead
    ggplot2::autoplot(z) # plot forecast

![](paper_files/figure-markdown_strict/1-ridge2-uni-1.png)

`EuStocksLogReturns` contains the daily log-returns of major European
stock indices, Germany DAX (Ibis), Switzerland SMI, France CAC, and UK
FTSE, with 1860 observations. Only the first 100 dates of DAX index are
used in the example below.

    data(EuStockMarkets)
    EuStocks <- ts(EuStockMarkets[1:100, ],
                   start = start(EuStockMarkets),
                   frequency = frequency(EuStockMarkets)) # original data
    EuStocksLogReturns <- ahead::getreturns(EuStocks, type = "log") # obtain log-returns
    res <- ahead::ridge2f(EuStocksLogReturns[, "DAX"], h = 20L,
                            show_progress = FALSE)
    ggplot2::autoplot(res) # plot forecast

![](paper_files/figure-markdown_strict/1-ridge2-uni-2-1.png)

## 2 - 2 - 2 Use `ahead::dynrmf` for univariate time series forecasting

`fdeaths` is used in this example too. The default model used by
`ahead::dynrmf` is an automated ridge regression (automatic choice of
the regularization parameter using Leave-One-Out cross-validation, see
Bergmeir, Hyndman, and Koo (2018)):

**- Forecasting with `randomForest::randomForest`**

    # Plotting forecasts
    # With a Random Forest regressor, an horizon of 20, 
    # and a 95% prediction interval
    fit_rf <- dynrmf(fdeaths, h=20, level=95, fit_func = randomForest::randomForest,
          fit_params = list(ntree = 50), predict_func = predict)
    ggplot2::autoplot(fit_rf)

![](paper_files/figure-markdown_strict/4-dynrmf-example-1.png)

Check in-sample residuals:

    forecast::checkresiduals(fit_rf)

![](paper_files/figure-markdown_strict/4-dynrmf-example-resids-1.png)

    ## 
    ##  Ljung-Box test
    ## 
    ## data:  Residuals from DynRM 1,1[12]
    ## Q* = 9.8649, df = 12, p-value = 0.6278
    ## 
    ## Model df: 0.   Total lags used: 12

**- Forecasting with `e1071::svm`** (Support Vector Machines)

    # With a Support Vector Machine regressor, an horizon of 20, 
    # and a 95% prediction interval
    fit_svm <- dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
    fit_params = list(kernel = "linear"), predict_func = predict)
    ggplot2::autoplot(fit_svm)

![](paper_files/figure-markdown_strict/5-dynrmf-example-1.png)

Check in-sample residuals:

    forecast::checkresiduals(fit_svm)

![](paper_files/figure-markdown_strict/5-dynrmf-example-resids-1.png)

    ## 
    ##  Ljung-Box test
    ## 
    ## data:  Residuals from DynRM 1,1[12]
    ## Q* = 27.351, df = 12, p-value = 0.006875
    ## 
    ## Model df: 0.   Total lags used: 12

**- Use of an external regressor** (trend)

`AirPassengers`’s been widely tested in the specialized literature,
because it has a trend, a seasonality, and is heteroskedastic
(non-constant variance).

    h <- 20L
    res6 <- ahead::dynrmf(AirPassengers, xreg_fit = 1:length(AirPassengers),
                           xreg_predict = (length(AirPassengers)+1):(length(AirPassengers)+h), 
                          h=h)
    autoplot(res6)

![](paper_files/figure-markdown_strict/unnamed-chunk-1-1.png)

## 2 - 2 - 2 Use `ahead::ridge2f` for multivariate time series forecasting

The `insurance` dataset available in Hyndman and Athanasopoulos (2013)
contains monthly quotations and monthly television advertising
expenditure for a US insurance company. January 2002 to April 2005. A
fast calibration of `ahead::ridge2f` uses a remarkable result available
for linear models’ Leave-One-Out Cross Validation (LOOCV):

$$
LOOCV error = \frac{1}{n} \sum\_{i=1}^n\left(y\_i-\hat{f}\_{-i}\left(\mathbf{z}\_i\right)\right)^2 = \frac{1}{n} \sum\_{i=1}^n\left(\frac{y\_i-\hat{f}\left(\mathbf{z}\_i\right)}{1-\mathbf{S}\_{i i}}\right)^2 
$$

Where *f̂*<sub>−*i*</sub>(**z**<sub>*i*</sub>) is a statistical learning
model *f*’s prediction without the *i*-th observation in the training
set. **z**<sub>*i*</sub> is the *i*-th observation, and
**S**<sub>*i**i*</sub> is the *i*-th diagonal element of the hat matrix
(a.k.a smoother, a matrix so that *ŷ* = **S***y*). **Keep in mind** that
[time series
cross-validation](https://thierrymoudiki.github.io/blog/2020/03/27/r/misc/crossval-2)
will give different results.

The LOOCV result can be extended further and approximated by
**Generalized Cross-Validation (GCV)** (Golub, Heath, and Wahba (1979)),
still with a closed-form formula available. GCV is used in
`ahead::ridge2f`. Indeed, even though `ridge2` is not – strictly
speaking – linear, it possesses the structure of a [ridge
regression](https://en.wikipedia.org/wiki/Ridge_regression) model with
[2 regularization parameters](https://www.mdpi.com/2227-9091/6/1/22); a
regularization parameter for the original explanative variables, and a
regularization parameter for new, engineered features. These engineered
features transform the linear model into a non-linear one.

In the following example, it’s worth mentioning that **only the 2
regularization parameters are calibrated**. Other model’s
hyperparameters such as the number of time series lags or the number of
nodes in the hidden layer are set to their default values (respectively
`1` and `5`).

    objective_function <- function(xx)
    {
        ahead::loocvridge2f(fpp::insurance,
                            h = 20,
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

    ## $par
    ## [1] -6.002722  1.247328
    ## 
    ## $value
    ## [1] 1.805636
    ## 
    ## $feval
    ## [1] 48
    ## 
    ## $restarts
    ## [1] 0
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## [1] "Successful convergence"

    print(proc.time()[3]-start)

    ## elapsed 
    ##  10.056

**Forecasting using the *optimal* regularization parameters**

    start <- proc.time()[3]
    res <- ahead::ridge2f(fpp::insurance, h = 20,
                          type_pi="blockbootstrap",
                          lambda_1=10^opt$par[1], # 'optimal' parameters
                          lambda_2=10^opt$par[2]) # 'optimal' parameters

    ##   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

    print(proc.time()[3]-start)

    ## elapsed 
    ##    0.25

    par(mfrow=c(2, 2))
    plot(res, "Quotes", type = "sims",
    main = "predictive simulations")
    plot(res, "TV.advert", type = "sims",
    main = "predictive simulations")
    plot(res, "Quotes", type = "dist",
    main = "prediction intervals")
    plot(res, "TV.advert", type = "dist",
    main = "prediction intervals")

![](paper_files/figure-markdown_strict/3-figure-2-real-world-1.png)

# References

Bergmeir, Christoph, Rob J Hyndman, and Bonsoo Koo. 2018. “A Note on the
Validity of Cross-Validation for Evaluating Autoregressive Time Series
Prediction.” *Computational Statistics & Data Analysis* 120: 70–83.

Brechmann, Eike Christian, and Ulf Schepsmeier. 2013. “Modeling
Dependence with c-and d-Vine Copulas: The r Package CDVine.” *Journal of
Statistical Software* 52: 1–27.

Efron, Bradley, and Robert Tibshirani. 1986. “Bootstrap Methods for
Standard Errors, Confidence Intervals, and Other Measures of Statistical
Accuracy.” *Statistical Science*, 54–75.

Golub, Gene H, Michael Heath, and Grace Wahba. 1979. “Generalized
Cross-Validation as a Method for Choosing a Good Ridge Parameter.”
*Technometrics* 21 (2): 215–23.

Goodfellow, Ian, Yoshua Bengio, and Aaron Courville. 2016. *Deep
Learning*. MIT press.

Hyndman, RJ, and G Athanasopoulos. 2013. “Forecasting: Principles and
Practice, OTexts. Org.” *URL: Https://Www. Otexts. Org/Fpp*.

Moudiki, Thierry, Frédéric Planchet, and Areski Cousin. 2018. “Multiple
Time Series Forecasting Using Quasi-Randomized Functional Link Neural
Networks.” *Risks* 6 (1): 22. <https://www.mdpi.com/2227-9091/6/1/22>.

Nagler, Thomas, Ulf Schepsmeier, Jakob Stoeber, Eike Christian
Brechmann, Benedikt Graeler, and Tobias Erhardt. 2023. *VineCopula:
Statistical Inference of Vine Copulas*.
<https://CRAN.R-project.org/package=VineCopula>.

Silverman, Bernard W. 2018. *Density Estimation for Statistics and Data
Analysis*. Routledge.

Vovk, Vladimir, Alexander Gammerman, and Glenn Shafer. 2005.
*Algorithmic Learning in a Random World*. Vol. 29. Springer.
