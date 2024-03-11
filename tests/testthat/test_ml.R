# context("Tests on `mlf`")
#
# # 1 - Code ----
#
# library(testthat)
# library(ahead)
#
# set.seed(123)
# x <- ts(rnorm(10))
# res <- ahead::embedc(x, 1)
#
# testthat::test_that("1 - tests on mlf", {
#   print("running tests on mlf")
#   expect_equal(res, embed(x, 1))
# })
