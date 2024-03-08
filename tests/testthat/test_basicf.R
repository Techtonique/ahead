library(testthat)
library(ahead)

set.seed(123)
x <- ts(matrix(rnorm(100), ncol=2))
z <- matrix(rnorm(100), ncol=2)

res1 <- ahead::basicf(x)
res2 <- ahead::basicf(x, type_pi = "bootstrap", B=10)
res3 <- ahead::basicf(x, type_pi = "blockbootstrap", B=10,
                      block_length = 4, show_progress = FALSE)
res4 <- ahead::basicf(z)
res5 <- ahead::basicf(x,
                      type_pi = "blockbootstrap",
                      B=10,
                      show_progress = FALSE)
res6 <- ahead::basicf(x,
                      type_pi = "movingblockbootstrap",
                      B=10,
                      show_progress = FALSE)


testthat::test_that("1 - tests on basicf", {
  print("running tests on basicf")
  expect_equal(as.numeric(round(res1$mean[1, 1], 2)), 0.03)
  expect_equal(as.numeric(round(res2$mean[1, 1], 2)), 0.03)
  expect_equal(as.numeric(round(res3$mean[1, 1], 2)), 0.03)
  expect_equal(as.numeric(round(res4$mean[1, 1], 2)), -0.25)
  expect_equal(as.numeric(round(res5$mean[1, 1], 2)), 0.03)
  expect_equal(as.numeric(round(res6$mean[1, 1], 2)), 0.03)
})
