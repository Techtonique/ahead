context("Tests on `dynrmf`")

# 1 - Code ----

library(testthat)
library(ahead)

set.seed(123)
x <- ts(rnorm(100))

res1 <- ahead::dynrmf(x)

testthat::test_that("1 - tests on dynrm", {
  expect_equal(as.numeric(round(res1$mean[1], 2)), 0.1)
})
