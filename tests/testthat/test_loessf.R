context("Tests on `loessf`")

# 1 - Code ----

library(testthat)
library(ahead)
require(e1071)

set.seed(123)
x <- ts(rnorm(50))

res1 <- ahead::loessf(x)
res2 <- ahead::loessf(x, type_aggregation = "mean")
res3 <- ahead::loessf(x, type_aggregation = "median")
res4 <- ahead::loessf(x, type_aggregation = "median",
                      type_pi = "bootstrap")
res5 <- ahead::loessf(x, type_aggregation = "median",
                      b = 3, type_pi = "blockbootstrap")
res6 <- ahead::loessf(x, type_aggregation = "median",
                      type_pi = "movingblockbootstrap")


testthat::test_that("1 - tests on dynrm", {
  print("running tests on dynrm")
  expect_equal(as.numeric(round(res1$mean[1], 2)), -0.24)
  expect_equal(as.numeric(round(res2$mean[1], 2)), -0.24)
  expect_equal(as.numeric(round(res3$mean[1], 2)), -0.35)
  expect_equal(as.numeric(round(res4$mean[1], 2)), -0.35)
  expect_equal(as.numeric(round(res5$mean[1], 2)),  -0.25)
  expect_equal(as.numeric(round(res6$mean[1], 2)), -0.42)
})
