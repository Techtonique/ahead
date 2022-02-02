context("Tests on `ridge2f`")


library(testthat)
library(ahead)

set.seed(123)
x <- ts(matrix(rnorm(100), ncol = 5))
res1 <- ahead::ridge2f(x)
res2 <- ahead::ridge2f(x, type_pi = "bootstrap", B=5)

# tests -----

test_that("tests on mean and lower", {
  expect_equal(as.numeric(round(res1$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res2$mean[1, 1], 2)), 1.01)
  expect_equal(as.numeric(round(res1$lower[1, 1], 2)), -0.34)
  expect_equal(as.numeric(round(res2$lower[1, 1], 2)), 0.31)
})
