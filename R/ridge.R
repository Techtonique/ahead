#' Fit Ridge regression
#'
#' @param x matrix of examples
#' @param y a vector, the response
#' @param lambda regularization parameters
#'
#' @return a list, an object of class 'ridge'
#'
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' n <- 100 ; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#'
#' fit_obj <- ahead::ridge(X, y)
#'
#' par(mfrow=c(1, 2))
#'
#' matplot(log(fit_obj$lambda), t(fit_obj$coef), type = 'l',
#' main="coefficients \n f(lambda)")
#'
#' plot(log(fit_obj$lambda), fit_obj$GCV, type='l',
#' main="GCV error")
#'
ridge <- function(x, y, lambda=10^seq(-10, 10,
                                          length.out = 100))
{
  # adapted from MASS::lm.ridge
  x <- as.matrix(x)
  y <- as.vector(y)
  nlambda <- length(lambda)

  ym <- mean(y)
  centered_y <- y - ym

  x_scaled <- base::scale(x)
  attrs <- attributes(x_scaled)
  X <- as.matrix(x_scaled[,])

  Xs <- La.svd(X)
  rhs <- crossprod(Xs$u, centered_y)
  d <- Xs$d
  nb_di <- length(d)
  div <- d ^ 2 + rep(lambda, rep(nb_di, nlambda))
  a <- drop(d * rhs) / div
  dim(a) <- c(nb_di, nlambda)
  n <- nrow(X)

  coef <- crossprod(Xs$vt, a)
  colnames(coef) <- lambda
  centered_y_hat <- X %*% coef

  fitted_values <- drop(ym +  centered_y_hat)
  if (length(lambda) > 1)
  {
    colnames(fitted_values) <- lambda
  }
  residuals <- centered_y - centered_y_hat
  GCV <- colSums(residuals^2)/(n - colSums(matrix(d^2/div, nb_di)))^2
  BIC <- n*log(colMeans(residuals^2)) + (ncol(X) + 2)*log(n)

  out <- list(
      coef = drop(coef),
      ym = ym,
      xm = attrs$`scaled:center`,
      xsd = attrs$`scaled:scale`,
      lambda = lambda,
      best_lam = lambda[which.min(GCV)],
      fitted_values = fitted_values,
      residuals = drop(centered_y - centered_y_hat),
      GCV = GCV,
      BIC = BIC,
      x = x,
      y = y
    )

  return(structure(out, class = "ridge"))
}


#' Predict from Ridge regression
#'
#' @param object object fitted by function `ridge`
#' @param newx new examples
#'
#' @return predicted values for \code{newx}
#' @export
#'
#' @examples
#'
#' set.seed(123)
#'
#' n <- 100 ; p <- 2
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#'
#' fit_obj <- ahead::ridge(X, y)
#'
#' n_test <- 10
#'
#' predict(fit_obj, newx=matrix(rnorm(n_test * p), n_test, p))
#'
#'
predict.ridge <- function(object, newx)
{
  if (length(object$lambda) > 1)
  {
    res <- try(drop(base::scale(newx, center=object$xm,
                                scale=object$xsd)%*%object$coef[,which.min(object$GCV)] + object$ym),
               silent = TRUE)
    if (inherits(res, "try-error"))
    {
      res <- try(drop(base::scale(newx, center=object$xm,
                                  scale=object$xsd)%*%object$coef[which.min(object$GCV)] + object$ym),
                 silent = TRUE)
      return(res)
    } else {
      return(res)
    }
  }  else {
    return(drop(base::scale(newx, center=object$xm,
                            scale=object$xsd)%*%object$coef + object$ym))
  }
}


#
# predict.ridge <- function(object, newx, cv=TRUE)
# {
#
#   if (cv){
#     res <- try(drop(base::scale(newx, center=object$xm,
#                                 scale=object$xsd)%*%object$coef[,which.min(object$GCV)] + object$ym),
#                silent = TRUE)
#     if (inherits(res, "try-error"))
#     {
#       res <- try(drop(base::scale(newx, center=object$xm,
#                                   scale=object$xsd)%*%object$coef[which.min(object$GCV)] + object$ym),
#                  silent = TRUE)
#       return(res)
#     } else {
#       return(res)
#     }
#   }
#
#   return(base::scale(newx, center=object$xm,
#                      scale=object$xsd)%*%object$coef + object$ym)
# }

