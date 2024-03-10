library(testthat)
library(ahead)

set.seed(123)
z <- matrix(rnorm(100), ncol=2)
x <- ts(z)

res1 <- ahead::varf(x)
res2 <- ahead::varf(z)

testthat::test_that("1 - tests on varf", {
  print("running tests on varf")
  expect_equal(as.numeric(round(res1$mean[1, 1], 2)), 0.04)
  expect_equal(as.numeric(round(res2$mean[1, 1], 2)), 0.04)
})
