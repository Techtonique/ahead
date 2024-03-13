context("Tests on `ridge2f`")

# 1 - Code ----

library(testthat)
library(ahead)

set.seed(123L)
x <- ts(rnorm(100L))

res1 <- ahead::armagarchf(x, h = 5,
                       level = 95,
                       B = 10L,
                       cl = 1L,
                       dist = "student",
                       seed = 1231)
res2 <- ahead::armagarchf(x, h = 6,
                       level = 80,
                       B = 10L,
                       cl = 1L,
                       dist = "gaussian",
                       seed = 1232)
res3 <- ahead::armagarchf(x, h = 7,
                       level = 95,
                       B = 10L,
                       cl = 1L,
                       dist = "student",
                       seed = 1233)
res4 <- ahead::armagarchf(x, h = 8,
                       level = 80,
                       B = 10L,
                       cl = 1L,
                       dist = "gaussian",
                       seed = 1234)
res5 <- ahead::armagarchf(x, h = 9,
                       level = 95,
                       B = 10L,
                       cl = 1L,
                       dist = "student",
                       seed = 1235)
res6 <- ahead::armagarchf(x, h = 5,
                          level = 95,
                          B = 10L,
                          cl = 2L,
                          dist = "student",
                          seed = 1231)

# 2 - tests -----

testthat::test_that("1 - tests on types of residuals' simulations", {
  expect_equal(as.numeric(round(res1$mean[5], 2)), -0.24)
  expect_equal(as.numeric(round(res2$mean[5], 2)), 0.95)
  expect_equal(as.numeric(round(res3$lower[5])), -1)
  expect_equal(as.numeric(round(res4$lower[5], 2)), -0.72)
  expect_equal(as.numeric(round(res5$mean[5], 2)), -0.18)
  expect_equal(as.numeric(round(res6$mean[5], 2)), -0.24)
})

