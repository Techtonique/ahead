
# y <- AirPassengers;
# obj <- model_fit(y, lags=2, obj = ahead::ridge);
# (res <- model_predict(obj, h=20L));
# plot(res, type = 'l')

# y <- Nile;
# obj <- model_fit(y, lags=2, obj = ahead::ridge);
# (res <- model_predict(obj, h=20L));
# plot(res, type = 'l')
#
# y <- USAccDeaths;
# obj <- model_fit(y, lags=2, obj = ahead::ridge);
# (res <- model_predict(obj, h=20L));
# plot(res, type = 'l')
#
# y <- mdeaths;
# obj <- model_fit(y, lags=2, obj = ahead::ridge);
# (res <- model_predict(obj, h=20L));
# plot(res, type = 'l')
#
# y <- austres
# obj <- model_fit(y, lags=2, obj = ahead::ridge);
# (res <- model_predict(obj, h=20L));
# plot(res, type = 'l')


# 0 - inputs and parameters -----
# For uncertainty:
# Gaussian
# Independant bootstrap
# MB bootstrap
# circular block bootstrap
# rgaussiandens
# rsurr
# rboot
model_fit <- function(y,
                      lags = 1L,
                      obj = randomForest::randomForest,
                      ...) {
  emb_y <- ahead::embedc(y, lags + 1L)
  p <- ncol(emb_y)
  if (lags == 1L) {
    x <- matrix(emb_y[,-p], ncol = 1)
  } else {
    x <- emb_y[,-p]
  }
  target <- ts(emb_y[, p], start=start(y),
               frequency=frequency(y))
  out <- list()
  out$lags <- lags
  out$x <- y
  out$y <- target
  out$reg <- x
  out$obj <- obj(out$reg, out$y, ...)
  return(out)
}

model_predict <- function(obj, h = 5, ...) {
  lags <- obj$lags
  p <- ncol(obj$reg)
  target <- obj$y
  tspy <- tsp(obj$x)
  start_preds <- tspy[2] + 1 / tspy[3]
  freq_preds <- frequency(obj$x)
  if (lags <= 1L) {
    newx <- matrix(target, ncol = 1)
  } else {
    newx <- cbind(obj$reg[, (p - lags + 2):p], target)
    colnames(newx) <- NULL
  }

  target <- predict(obj$obj, newx, ...)
  for (i in 2:h) {
    if (lags <= 1L) {
      newx <- matrix(target, ncol = 1)
    } else {
      newx <- cbind(newx[, (p - lags + 2):p], target)
      colnames(newx) <- NULL
    }
    target <- predict(obj$obj, newx, ...)
  }
  return(target)
}


# y <- AirPassengers; print(y); obj <- model_fit(y, lags=2); (res <- model_predict(obj, h=20L))
