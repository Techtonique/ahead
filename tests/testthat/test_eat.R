context("Tests on `eatf`")

# 1 - Code ----

library(testthat)
library(ahead)

set.seed(123)
x <- ts(rnorm(100))

res1 <- ahead::eatf(x)

testthat::test_that("1 - tests on eat", {
  expect_equal(as.numeric(round(res1$mean[1], 2)), 0.1)
})
