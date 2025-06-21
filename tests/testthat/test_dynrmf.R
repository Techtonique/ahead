context("Tests on `dynrmf`")

# 1 - Code ----

library(testthat)
library(ahead)
require(e1071)

set.seed(123)
x <- ts(rnorm(100))

res5 <- ahead::dynrmf(x, fit_func = e1071::svm,
       fit_params = list(kernel = "linear"),
       predict_func = predict)
res7 <- ahead::dynrmf(x, xreg_fit = 1:length(x),
                      xreg_predict = (length(x)+1):(length(x)+5))
x <- AirPassengers
res6 <- ahead::dynrmf(x, xreg_fit = 1:length(x),
                      xreg_predict = (length(x)+1):(length(x)+5))
res8 <- ahead::dynrmf(abs(x), lambda="auto")

testthat::test_that("1 - tests on dynrm", {
  print("running tests on dynrm")
  expect_equal(as.numeric(round(res5$mean[1], 2)), 0.08)
  expect_equal(as.numeric(round(res6$mean[1], 2)), 455.82)
  expect_equal(as.numeric(round(res7$mean[1], 2)), 0.1)
  expect_equal(as.numeric(round(res8$mean[1], 2)), 454.44)
})
