
# ahead

Univariate and multivariate time series forecasting, with uncertainty quantification.

[![Documentation](https://img.shields.io/badge/documentation-is_here-green)](https://techtonique.github.io/ahead/index.html)
![Downloads](https://r-packages.techtonique.net/badges/downloads/last-month/ahead.svg)
![Total Downloads](https://r-packages.techtonique.net/badges/downloads/grand-total/ahead.svg?color=brightgreen)


### Installation

- __1st method__, released version: from [Techtonique's R package's repository](https://r-packages.techtonique.net)

    In R console:
    
    ```R
    options(repos = c(
                    techtonique = "https://r-packages.techtonique.net",
                    CRAN = "https://cloud.r-project.org"
                ))            
        
    install.packages("ahead")
    ```

- __2nd method__: from Github

    In R console:
    
    ```R
    devtools::install_github("Techtonique/ahead")
    ```

### Demo

For univariate and multivariate time series.

#### 1 - Univariate time series 

##### 1 - 1 Example 1: with `dynrmf` (type `?dynrmf` in R console for more details) and Random Forest

```R
 require(fpp2)
 
 par(mfrow=c(3, 2))
 plot(dynrmf(USAccDeaths, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict))
 plot(dynrmf(AirPassengers, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict))
 plot(dynrmf(lynx, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict))
 plot(dynrmf(WWWusage, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict))
 plot(dynrmf(Nile, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict))
 plot(dynrmf(fdeaths, h=20, level=95, fit_func = randomForest::randomForest,
      fit_params = list(ntree = 50), predict_func = predict))
```      


##### 1 - 2 Example 2:  with `dynrmf` and Support Vector Machines

```R
 require(e1071)
 
 par(mfrow=c(2, 2))
 plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
 fit_params = list(kernel = "linear"), predict_func = predict))
 plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
 fit_params = list(kernel = "polynomial"), predict_func = predict))
 plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
 fit_params = list(kernel = "radial"), predict_func = predict))
 plot(dynrmf(fdeaths, h=20, level=95, fit_func = e1071::svm,
 fit_params = list(kernel = "sigmoid"), predict_func = predict))
```

For more examples on `dynrmf`, you can read this  [blog post](https://thierrymoudiki.github.io/blog/2021/12/20/r/forecasting/ahead-more-examples).

#### 2 - Multivariate time series

With `ridge2f` (type `?ridge2f` in R console for more details), the model from : 

 Moudiki, T., Planchet, F., & Cousin, A. (2018).
 Multiple time series forecasting using quasi-randomized
 functional link neural networks. Risks, 6(1), 22.

```R
 require(fpp2)

 print(ahead::ridge2f(fpp2::insurance)$mean)

 res <- ahead::ridge2f(fpp2::insurance, lags=2)
 par(mfrow=c(1, 2))
 plot(res, "Quotes")
 plot(res, "TV.advert")
```

### Contributing

Your contributions are welcome. Please, make sure to __read__ the [Code of Conduct](CONTRIBUTING.md) first.

### Note to self

```bash
git remote set-url origin https://username:your_generated_token@github.com/xxx/repo.git
```

```bash
git remote set-url origin https://MY_GITHUB_TOKEN@github.com/Techtonique/ahead.git
```

### License

[BSD 3-Clause](LICENSE) © Thierry Moudiki, 2019. 

