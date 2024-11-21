# version 0.14.0

- Add `fit_func` and `predict_func` for custom fitting and prediction functions of `ahead::dynrmf` (using `caret` Machine Learning).
- Add forecasting combinations based on ForecastComb, adding Ridge and Elastic Net to the mix.

# version 0.11.0

- Include tests (90% coverage). After cloning, run: 

```R
install.packages("covr")
covr::report()
```

# version 0.10.0

- Univariate forecasting for `ridge2f`. 
See https://thierrymoudiki.github.io/blog/2024/02/26/python/r/julia/ahead-v0100.
- Fast calibration for `ridge2f` (univariate and multivariate case). 
See https://thierrymoudiki.github.io/blog/2024/02/26/python/r/julia/ahead-v0100.

# version 0.9.0

- progress bars for bootstrap (independent, circular block, moving block)

# version 0.8.0

- empirical marginals for R-Vine copula simulation 
- risk-neutralize simulations

# version 0.7.0

- moving block bootstrap in `ridge2f`, `basicf` and `loessf`, in addition to circular block bootstrap from 0.6.2
- adjust R-Vine copulas on residuals for `ridge2f` simulation
- new plots for simulations see (new) vignettes
- split conformal prediction intervals (**very very experimental** and basic right now, too conservative)
- `Depends` and selective `Imports` (beneficial to Python and rpy2 for installation time?)
- `getsimulations` extracts simulations from a given time series (from `ridge2f` and `basicf`)
- `getreturns` extracts returns/log-returns from multivariate time series
- `splitts` splits time series using a proportion of data

# version 0.6.2

- Add Block Bootstrap to `ridge2f`
- Add external regressors to `ridge2f`
- Add clustering to `ridge2f`
- Add Block Bootstrap to `loessf`
- Create new vignettes for `ridge2f` and `loessf`

# version 0.6.1

- Align version with Python's 
- Temporarily remove dependency with `cclust`

# version 0.6.0

- Include basic methods: mean forecast, median forecast, random walk forecast

# version 0.5.0

- add dropout regularization to `ridge2f`
- parallel execution for `type_pi == bootstrap` in `ridge2f` (done in R /!\, experimental)
- preallocate matrices for `type_forecast == recursive` in `ridge2f`


# version 0.4.2

- new attributes mean, lower bound, upper bound forecast as numpy arrays


# version 0.4.1

- use `get_frequency` to get series frequency as a number
- create a function `get_tscv_indices` for getting time series cross-validation indices
