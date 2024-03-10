context("Tests on `ridge2f`")

# 1 - Code ----

library(testthat)
library(ahead)

set.seed(123L)
z <- matrix(rnorm(100L), ncol = 5L)
x <- ts(z)

# 1 - 1 type of residuals 'simulation' ----

res1 <- ahead::ridge2f(x, type_pi = "gaussian", B = 5L)
res37 <- ahead::ridge2f(z, type_pi = "gaussian", B = 5L)
res2 <- ahead::ridge2f(x, type_pi = "bootstrap", B = 5L)
res3 <- ahead::ridge2f(x, type_pi = "bootstrap", B = 5L)
res4 <- ahead::ridge2f(x, type_pi = "blockbootstrap", B = 5L)
res5 <- ahead::ridge2f(x, type_pi = "movingblockbootstrap", B = 5L)
res6 <- ahead::ridge2f(x, type_pi = "rvinecopula", B = 5L)
res41 <- ahead::ridge2f(x, type_pi = "rvinecopula", B = 5L,
                        margins = "empirical")
res43 <- ahead::ridge2f(x, type_pi = "rvinecopula", B = 5L,
                        margins = "empirical", cl=2L)
res24 <- ahead::ridge2f(x,
                        type_pi = "bootstrap",
                        B = 5L,
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

res16 <- ahead::ridge2f(x, centers = 2L,
                        type_clustering = "kmeans")
res17 <- ahead::ridge2f(x, centers = 2L,
                        type_clustering = "hclust")

# 1 - 5 xreg ----

xreg <- 1:nrow(x)
xreg2 <- cbind(1:nrow(x), rnorm(nrow(x)))
xreg3 <- cbind(1:nrow(x), rnorm(nrow(x)))
colnames(xreg3) <- c("col1", "col2")
res18 <- ahead::ridge2f(x, xreg = xreg)
res19 <- ahead::ridge2f(x, xreg = xreg2)
res20 <- ahead::ridge2f(x, xreg = xreg3)

# 1 - 6 xreg with clustering ----

res21 <-
  ahead::ridge2f(x,
                 xreg = xreg,
                 centers = 2L,
                 type_clustering = "kmeans")
res22 <-
  ahead::ridge2f(x,
                 xreg = xreg2,
                 centers = 2L,
                 type_clustering = "hclust")
res23 <-
  ahead::ridge2f(x,
                 xreg = xreg3,
                 centers = 2L,
                 type_clustering = "kmeans")

res21_ <-
  ahead::ridge2f(x,
                 xreg = xreg,
                 centers = 2L,
                 type_clustering = "kmeans",
                 show_progress = TRUE,
                 type_pi = "bootstrap",
                 cl = 2L)
res22_ <-
  ahead::ridge2f(x,
                 xreg = xreg2,
                 centers = 2L,
                 type_clustering = "hclust",
                 show_progress = TRUE,
                 type_pi = "bootstrap",
                 cl = 2L)
res23_ <-
  ahead::ridge2f(x,
                 xreg = xreg3,
                 centers = 2L,
                 type_clustering = "kmeans",
                 show_progress = TRUE,
                 type_pi = "bootstrap",
                 cl = 2L)

# 1 - 7 direct forecasting ----

res25 <- ahead::ridge2f(x, type_forecast = "direct")
res26 <- ahead::ridge2f(
  x,
  type_forecast = "direct",
  type_pi = "bootstrap",
  B = 5L,
  show_progress = FALSE
)
res26_ <- ahead::getsimulations(res26, 1)
res26__ <- ahead::getsimulations(res26, 1, transpose = TRUE)

# 1 - 8 univariate with xreg ----

x <- fdeaths
xreg <- ahead::createtrendseason(x) # add seasonality and trend
res27 <-
  ahead::ridge2f(x, xreg = xreg, h = 5L) # forecasting h-steps ahead


# 1 - 9 dropout ----

set.seed(123L)
x <- ts(matrix(rnorm(100L), ncol = 5L))
res28 <- ahead::ridge2f(x, dropout = 0.1)


# 1 - 10 get returns ----
res29 <- getreturns(EuStockMarkets)
res30 <- getreturns(EuStockMarkets,
                          type = "log")

# 1 - 11 loocv ----

set.seed(123L)
x <- ts(matrix(rnorm(50L), ncol = 2L))
objective_function <- function(xx)
{
  ahead::loocvridge2f(x,
                      lambda_1=10^xx[1],
                      lambda_2=10^xx[2],
                      show_progress = FALSE,
  )$loocv
}
(opt <- dfoptim::nmkb(fn=objective_function,
                      lower=c(-10,-10),
                      upper=c(10,10),
                      par=c(0.1, 0.1)))
res31 <- ahead::ridge2f(x,
                      lambda_1=10^opt$par[1], # 'optimal' parameters
                      lambda_2=10^opt$par[2]) # 'optimal' parameters


# 1 - 12 plot ----

# moving block bootstrap
res32 <- ahead::ridge2f(x, type_pi = "bootstrap",
                        B=5, show_progress = FALSE)

res33 <- plot(res32,
              selected_series="Series 1",
              type = "pi")

res34 <- plot(res32,
              selected_series="Series 1",
              type = "dist")

res35 <- plot(res32,
              selected_series="Series 1",
              type = "sims")

# 1 - 13 parallel exec ----

set.seed(123L); x <- ts(matrix(rnorm(50L), ncol = 2L))
res36 <- ahead::ridge2f(x, type_pi = "movingblockbootstrap",
                        B = 10L, cl=2L)
res38 <- ahead::ridge2f(x, type_pi = "movingblockbootstrap",
                        B = 10L, cl=2L,
                        show_progress = FALSE)

# 1 - 14 ym ----

data(EuStockMarkets)
EuStocks <- ts(EuStockMarkets[1:100, ],
               start = start(EuStockMarkets),
               frequency = frequency(EuStockMarkets))
EuStocksLogReturns <- ahead::getreturns(EuStocks, type = "log")
ym <- c(0.03013425, 0.03026776, 0.03040053, 0.03053258,
        0.03066390, 0.03079450, 0.03092437)
freq <- frequency(EuStocksLogReturns)
(start_preds <- tsp(EuStocksLogReturns)[2] + 1 / freq)
(ym <- stats::ts(ym,
                 start = start_preds,
                 frequency = frequency(EuStocksLogReturns)))
res39 <- ahead::ridge2f(EuStocksLogReturns, h = 7L,
                        type_pi = 'bootstrap',
                        B = 10L, ym = ym,
                        show_progress = FALSE)

# 2 - tests -----

# 2 - 1 type of residuals 'simulation' ----

testthat::test_that("1 - tests on types of residuals' simulations", {
  print("running tests on types of residuals' simulations")
  expect_equal(as.numeric(round(res1$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res37$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res2$mean[1, 1], 2)), 1.01)
  expect_equal(as.numeric(round(res1$lower[1, 1], 2)),-0.34)
  expect_equal(as.numeric(round(res2$lower[1, 1], 2)), 0.31)
  expect_equal(as.numeric(round(res3$mean[1, 1], 2)), 1.01)
  expect_equal(as.numeric(round(res4$mean[1, 1], 2)), 1.46)
  expect_equal(as.numeric(round(res3$lower[1, 1], 2)), 0.31)
  expect_equal(as.numeric(round(res4$lower[1, 1], 2)), 0.55)
  expect_equal(as.numeric(round(res41$mean[1, 1], 2)), 1.34)
  expect_equal(as.numeric(round(res43$mean[1, 1], 2)), 1.34)
  expect_equal(as.numeric(round(res24$mean[1, 1], 2)), 1.01)
  expect_equal(as.numeric(round(res27$mean[1], 2)), 727.05)
  expect_equal(as.numeric(round(res28$mean[1], 2)), 0.8)
})

# 2 - 2 type of residuals nodes' simulation' ----

testthat::test_that("2 - tests on types of nodes' simulations", {
  print("running tests on types of nodes' simulations")
  expect_equal(as.numeric(round(res7$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res8$mean[1, 1], 2)), 0.32)
  expect_equal(as.numeric(round(res9$mean[1, 1], 2)), 0.83)
})

# 2 - 3 type of activation function ----

testthat::test_that("3 - tests on types of activation function", {
  print("running tests on types of activation function")
  expect_equal(as.numeric(round(res10$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res11$mean[1, 1], 2)), 0.19)
  expect_equal(as.numeric(round(res12$mean[1, 1], 2)), 0.22)
  expect_equal(as.numeric(round(res13$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res14$mean[1, 1], 2)), 1.16)
  expect_equal(as.numeric(round(res15$mean[1, 1], 2)), 0.32)

})

# 2 - 4 clustering ----

testthat::test_that("4 - tests on clustering", {
  print("running tests on clustering")
  expect_equal(as.numeric(round(res16$mean[1, 1], 2)),-0.43)
  expect_equal(as.numeric(round(res17$mean[1, 1], 2)),-0.15)
})

# 2 - 5 xreg ----

testthat::test_that("5 - tests on xreg", {
  print("running tests on xreg")
  expect_equal(as.numeric(round(res18$mean[1, 1], 2)),-0.03)
  expect_equal(as.numeric(round(res19$mean[1, 1], 2)),-0.16)
  expect_equal(as.numeric(round(res20$mean[1, 1], 2)),-0.71)
})

# 2 - 6 xreg with clustering ----

testthat::test_that("6 - tests on xreg and clustering", {
  print("running tests on xreg and clustering")
  expect_equal(as.numeric(round(res21$mean[1, 1], 2)),-0.13)
  expect_equal(as.numeric(round(res22$mean[1, 1], 2)), 0.04)
  expect_equal(as.numeric(round(res23$mean[1, 1], 2)),-0.76)
  expect_equal(as.numeric(round(res21_$mean[1, 1], 2)),-0.26)
  expect_equal(as.numeric(round(res22_$mean[1, 1], 2)), -0.04)
  expect_equal(as.numeric(round(res23_$mean[1, 1], 2)), -0.83)
})


# 2 - 7 direct forecasting ----

testthat::test_that("7 - tests on direct forecasting", {
  print("running tests on direct forecasting")
  expect_equal(as.numeric(round(res25$mean[1, 3], 2)), 0.74)
  expect_equal(as.numeric(round(res25$upper[1, 3], 2)), 2.32)
  expect_equal(as.numeric(round(res26$mean[1, 3], 2)), 0.86)
  expect_equal(as.numeric(round(res26$upper[1, 3], 2)), 1.41)
  expect_equal(as.numeric(round(res26_$series[1, 1], 2)), 0.45)
  expect_equal(as.numeric(round(res26__$series[1, 1], 2)), 0.45)
})


# 2 - 8 get returns and log-returns -----

testthat::test_that("8 - tests on returns", {
  print("running tests on returns")
  expect_equal(as.numeric(round(getreturns(EuStockMarkets)[1,"DAX"], 2)), -0.01)
  expect_equal(as.numeric(round(getreturns(EuStockMarkets, type="log")[1,"DAX"], 2)), -0.01)
})


# 2 - 9 loocv -----

testthat::test_that("9 - tests on loocv ridge2", {
  print("running tests on loocv ridge2")
  expect_equal(as.numeric(round(res31$mean[1, 1], 2)), -0.01)
})


# 2 - 10 plot -----

testthat::test_that("10 - tests on plots", {
  print("running tests on plots")
  expect_equal(res33, NULL)
  expect_equal(res34, NULL)
  expect_equal(res35, NULL)
})


# 2 - 11 parallel -----

testthat::test_that("11 - parallel exec", {
  print("running tests on parallel exec")
  expect_equal(as.numeric(round(res36$mean[1, 1], 2)),
               0.35)
  expect_equal(as.numeric(round(res38$mean[1, 1], 2)),
               0.35)
})

# 2 - 12 ym -----

testthat::test_that("11 - ym", {
  print("running tests on ym")
  expect_equal(as.numeric(round(res39$mean[1, 1], 4)),
               0.0026)
})
