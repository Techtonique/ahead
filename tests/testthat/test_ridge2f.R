context("Tests on `ridge2f`")

# 1 - Code ----

library(testthat)
library(ahead)

set.seed(123)
x <- ts(matrix(rnorm(100), ncol = 5))

# 1 - 1 type of residuals 'simulation' ----

res1 <- ahead::ridge2f(x, type_pi = "gaussian", B=5)
res2 <- ahead::ridge2f(x, type_pi = "bootstrap", B=5)
res3 <- ahead::ridge2f(x, type_pi = "bootstrap", B=5)
res4 <- ahead::ridge2f(x, type_pi = "blockbootstrap", B=5)
res5 <- ahead::ridge2f(x, type_pi = "movingblockbootstrap", B=5)
res6 <- ahead::ridge2f(x, type_pi = "rvinecopula", B=5)
res24 <- ahead::ridge2f(x, type_pi = "bootstrap", B=5,
                        show_progress = FALSE)

# 1 - 2 type of residuals nodes' simulation' ----

res7 <- ahead::ridge2f(x, nodes_sim = "sobol")
res8 <- ahead::ridge2f(x, nodes_sim = "halton")
res9 <- ahead::ridge2f(x, nodes_sim = "unif")

# 1 - 3 type of activation function ----

res10 <- ahead::ridge2f(x, activ = "relu")
res11 <- ahead::ridge2f(x, activ = "sigmoid")
res12 <- ahead::ridge2f(x, activ = "tanh")
res13 <- ahead::ridge2f(x, activ = "leakyrelu")
res14 <- ahead::ridge2f(x, activ = "elu")
res15 <- ahead::ridge2f(x, activ = "linear")

# 1 - 4 clustering ----

res16 <- ahead::ridge2f(x, centers = 2L, type_clustering = "kmeans")
res17 <- ahead::ridge2f(x, centers = 2L, type_clustering = "hclust")

# 1 - 5 xreg ----

xreg <- 1:nrow(x)
xreg2 <- cbind(1:nrow(x), rnorm(nrow(x)))
xreg3 <- cbind(1:nrow(x), rnorm(nrow(x)))
colnames(xreg3) <- c("col1", "col2")
res18 <- ahead::ridge2f(x, xreg = xreg)
res19 <- ahead::ridge2f(x, xreg = xreg2)
res20 <- ahead::ridge2f(x, xreg = xreg3)

# 1 - 6 xreg with clustering ----

res21 <- ahead::ridge2f(x, xreg = xreg, centers = 2L, type_clustering = "kmeans")
res22 <- ahead::ridge2f(x, xreg = xreg2, centers = 2L, type_clustering = "hclust")
res23 <- ahead::ridge2f(x, xreg = xreg3, centers = 2L, type_clustering = "kmeans")


# 2 - tests -----

# 2 - 1 type of residuals 'simulation' ----

testthat::test_that("1 - tests on types of residuals' simulations", {
  expect_equal(as.numeric(round(res1$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res2$mean[1, 1], 2)), 1.01)
  expect_equal(as.numeric(round(res1$lower[1, 1], 2)), -0.34)
  expect_equal(as.numeric(round(res2$lower[1, 1], 2)), 0.31)
  expect_equal(as.numeric(round(res3$mean[1, 1], 2)), 1.01)
  expect_equal(as.numeric(round(res4$mean[1, 1], 2)), 1.46)
  expect_equal(as.numeric(round(res3$lower[1, 1], 2)), 0.31)
  expect_equal(as.numeric(round(res4$lower[1, 1], 2)), 0.55)
  expect_equal(as.numeric(round(res24$mean[1, 1], 2)), 1.01)
})

# 2 - 2 type of residuals nodes' simulation' ----

testthat::test_that("2 - tests on types of nodes' simulations", {
  expect_equal(as.numeric(round(res7$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res8$mean[1, 1], 2)), 0.32)
  expect_equal(as.numeric(round(res9$mean[1, 1], 2)), 0.83)
})

# 2 - 3 type of activation function ----

testthat::test_that("3 - tests on types of activation function", {
  expect_equal(as.numeric(round(res10$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res11$mean[1, 1], 2)), 0.19)
  expect_equal(as.numeric(round(res12$mean[1, 1], 2)), 0.22)
  expect_equal(as.numeric(round(res13$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res14$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res15$mean[1, 1], 2)), 0.32)

})

# 2 - 4 clustering ----

testthat::test_that("4 - tests on clustering", {
  expect_equal(as.numeric(round(res16$mean[1, 1], 2)), -0.43)
  expect_equal(as.numeric(round(res17$mean[1, 1], 2)), -0.15)
})

# 2 - 5 xreg ----

testthat::test_that("4 - tests on xreg", {
  expect_equal(as.numeric(round(res18$mean[1, 1], 2)), -0.03)
  expect_equal(as.numeric(round(res19$mean[1, 1], 2)), -0.16)
  expect_equal(as.numeric(round(res20$mean[1, 1], 2)), -0.71)
})

# 2 - 5 xreg with clustering ----

testthat::test_that("5 - tests on xreg and clustering", {
  expect_equal(as.numeric(round(res21$mean[1, 1], 2)), -0.13)
  expect_equal(as.numeric(round(res22$mean[1, 1], 2)), 0.04)
  expect_equal(as.numeric(round(res23$mean[1, 1], 2)), -0.76)
})
