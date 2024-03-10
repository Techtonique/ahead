context("Tests on `eatf`")

# 1 - Code ----

library(testthat)
library(ahead)

set.seed(123)
x <- ts(rnorm(100))

# method = c("EAT", "E", "A", "T"),
# weights = rep(1/3, 3),
# type_pi = c("gaussian", "E", "A", "T")
res1 <- ahead::eatf(x)
res2 <- ahead::eatf(x, method = "E")
res3 <- ahead::eatf(x, method = "A")
res4 <- ahead::eatf(x, method = "T")
res5 <- ahead::eatf(x, method = "EAT")
res6 <- ahead::eatf(x, weights = c(0, 0.5, 0.5))
res7 <- ahead::eatf(x, weights = c(0.5, 0, 0.5))
res8 <- ahead::eatf(x, weights = c(0.5, 0.5, 0))
res9 <- ahead::eatf(x, type_pi = "E")
res10 <- ahead::eatf(x, type_pi = "A")
res11 <- ahead::eatf(x, type_pi = "T")

testthat::test_that("1 - tests on eat", {
  print("running tests on eat")
  expect_equal(as.numeric(round(res1$mean[1], 2)), 0.1)
  expect_equal(as.numeric(round(res2$mean[1], 2)), 0.1)
  expect_equal(as.numeric(round(res3$mean[1], 2)), 0)
  expect_equal(as.numeric(round(res4$mean[1], 2)), 0.22)
  expect_equal(as.numeric(round(res5$mean[1], 2)), 0.1)
  expect_equal(as.numeric(round(res6$mean[1], 2)), 0.11)
  expect_equal(as.numeric(round(res7$mean[1], 2)), 0.16)
  expect_equal(as.numeric(round(res8$mean[1], 2)), 0.05)
  expect_equal(as.numeric(round(res9$mean[1], 2)), 0.14)
  expect_equal(as.numeric(round(res10$mean[1], 2)), 0.1)
  expect_equal(as.numeric(round(res11$mean[1], 2)), 0.26)
})
