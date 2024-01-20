library(testthat)
library(ahead)

set.seed(123)
x <- ts(matrix(rnorm(100), ncol=2))

res1 <- ahead::varf(x)

testthat::test_that("1 - tests on dynrm", {
  expect_equal(as.numeric(round(res1$mean[1, 1], 2)), 0.04)
})
