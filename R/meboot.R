#' Maximum Entropy Bootstrap for Time Series using Rcpp
#'
#' Generates bootstrap replicates of a time series using the maximum entropy 
#' bootstrap algorithm with Rcpp implementation for improved performance.
#' This method is particularly useful for non-stationary time series and 
#' preserves the dependence structure of the original data.
#'
#' @param x A numeric vector or time series object to be bootstrapped.
#' @param reps Number of bootstrap replicates to generate (default: 999).
#' @param trim Controls tail behavior. Can be a single numeric value specifying 
#'   the trim proportion (e.g., 0.10 for 10% trimming), or a list with components:
#'   \describe{
#'     \item{trim}{Trim proportion for tail calculation (default: 0.10)}
#'     \item{xmin}{Lower bound for generated values (optional)}
#'     \item{xmax}{Upper bound for generated values (optional)}
#'   }
#' @param reachbnd Logical indicating whether to allow generated values to reach 
#'   the boundaries xmin and xmax (default: TRUE).
#' @param expand.sd Logical indicating whether to expand the standard deviation 
#'   of the ensemble (default: TRUE).
#' @param force.clt Logical indicating whether to force the central limit theorem 
#'   compliance by centering each replicate (default: TRUE).
#' @param scl.adjustment Logical indicating whether to adjust the scale of the 
#'   ensemble to match the original data's variance (default: FALSE).
#' @param sym Logical indicating whether to force symmetry in the maximum entropy 
#'   density (default: FALSE).
#' @param colsubj Deprecated parameter from original meboot (included for compatibility).
#' @param coldata Deprecated parameter from original meboot (included for compatibility).
#' @param coltimes Deprecated parameter from original meboot (included for compatibility).
#' @param ... Additional arguments passed to expansion functions.
#'
#' @return A list with the following components:
#'   \item{x}{Original time series data}
#'   \item{ensemble}{Matrix of bootstrap replicates (n x reps)}
#'   \item{xx}{Sorted original data}
#'   \item{z}{Intermediate points between sorted values}
#'   \item{dv}{Absolute differences between consecutive observations}
#'   \item{dvtrim}{Trimmed mean of differences}
#'   \item{xmin}{Lower bound used for generation}
#'   \item{xmax}{Upper bound used for generation}
#'   \item{desintxb}{Interval means satisfying mean-preserving constraint}
#'   \item{ordxx}{Ordering index of original data}
#'   \item{kappa}{Scale adjustment factor (if scl.adjustment = TRUE)}
#'
#' @details
#' The maximum entropy bootstrap algorithm generates replicates that:
#' \itemize{
#'   \item Preserve the dependence structure of the original time series
#'   \item Can handle non-stationary time series
#'   \item Satisfy the ergodic theorem and central limit theorem
#'   \item Maintain the mean and autocorrelation structure
#' }
#' 
#' The Rcpp implementation provides significant performance improvements 
#' over the original R implementation, especially for large datasets 
#' and many replications.
#'
#' @examples
#' \donttest{
#' # Basic usage with a time series
#' set.seed(123)
#' x <- ts(rnorm(100), start = c(2000, 1), frequency = 12)
#' result <- meboot(x, reps = 1000)
#' 
#' # Plot first few replicates
#' matplot(result$ensemble[, 1:5], type = "l", lty = 1)
#' lines(result$x, col = "black", lwd = 2)
#' 
#' # With custom bounds
#' result_bounded <- meboot(x, reps = 100, 
#'                              trim = list(trim = 0.1, xmin = -3, xmax = 3))
#' 
#' # With scale adjustment
#' result_scaled <- meboot(x, reps = 100, scl.adjustment = TRUE)
#' }
#'
#' @references
#' 
#' Vinod, H. D., & Lopez-de-Lacalle, J. (2009). Maximum entropy 
#' bootstrap for time series: The meboot R package. 
#' \emph{Journal of Statistical Software}, 29(5), 1-19.
#'
#' @export
#'
#' @keywords ts bootstrap
#' 
meboot <- function(x, reps = 999, 
                   trim = list(trim = 0.10, xmin = NULL, xmax = NULL), 
                   reachbnd = TRUE, expand.sd = TRUE, force.clt = TRUE,
                   scl.adjustment = FALSE, sym = FALSE, 
                   colsubj, coldata, coltimes, ...) {
  
  # Load Rcpp
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp package is required")
  }
  
  # Handle pdata.frame case (simplified)
  if ("pdata.frame" %in% class(x)) {
    warning("pdata.frame support not fully implemented in Rcpp version")
  }
  
  if (reps == 1 && isTRUE(force.clt)) {
    force.clt <- FALSE
    warning("force.clt was set to FALSE since the ensemble contains only one replicate.")
  }
  
  if (!is.list(trim)) {
    trimval <- trim
  } else {
    trimval <- if (is.null(trim$trim)) 0.1 else trim$trim
  }
  
  ptm1 <- proc.time()
  n <- length(x)
  
  # Sort the original data
  xx <- sort(x)
  ordxx <- order(x)
  
  # symmetry
  if (sym) {
    xxr <- rev(xx)
    xx.sym <- mean(xx) + 0.5 * (xx - xxr)
    xx <- xx.sym
  }
  
  # Compute intermediate points
  z <- (xx[-1] + xx[-n]) / 2
  
  # Compute limits
  dv <- abs(diff(as.numeric(x)))
  dvtrim <- mean(dv, trim = trimval)
  
  if (is.list(trim)) {
    xmin <- if (is.null(trim$xmin)) xx[1] - dvtrim else trim$xmin
    xmax <- if (is.null(trim$xmax)) xx[n] + dvtrim else trim$xmax
    
    if (!is.null(trim$xmin) || !is.null(trim$xmax)) {
      if (isTRUE(force.clt)) {
        expand.sd <- FALSE
        force.clt <- FALSE
        warning("expand.sd and force.clt were set to FALSE to enforce limits xmin/xmax.")
      }
    }
  } else {
    xmin <- xx[1] - dvtrim
    xmax <- xx[n] + dvtrim
  }
  
  # Warnings for limits
  if (is.list(trim)) {
    if (!is.null(trim$xmin) && trim$xmin > min(x))
      warning("trim$xmin may not be satisfied in replicates")
    if (!is.null(trim$xmax) && trim$xmax < max(x))
      warning("trim$xmax may not be satisfied in replicates")
  }
  
  # Compute interval means
  aux <- colSums(t(cbind(xx[-c(1,2)], xx[-c(1,n)], xx[-c((n-1),n)])) * c(0.25, 0.5, 0.25))
  desintxb <- c(0.75 * xx[1] + 0.25 * xx[2], aux, 0.25 * xx[n-1] + 0.75 * xx[n])
  
  # Generate ensemble using Rcpp
  # NumericVector p, int n, NumericVector z, NumericVector desintxb
  ensemble <- meboot_part_rcpp(x, reps, z, xmin, xmax, desintxb, reachbnd)
  
  # Apply ordering
  qseq <- apply(ensemble, 2, sort)
  ensemble[ordxx, ] <- qseq
  
  # Expand SD if requested
  if (expand.sd) {
    ensemble <- expand_sd_rcpp(x, ensemble, ...)
  }
  
  # Force CLT if requested (simplified implementation)
  if (force.clt) {
    x_mean <- mean(x)
    for (j in 1:reps) {
      ens_mean <- mean(ensemble[, j])
      ensemble[, j] <- ensemble[, j] - ens_mean + x_mean
    }
  }
  
  # Scale adjustment
  if (scl.adjustment) {
    zz <- c(xmin, z, xmax)
    v <- diff(zz^2) / 12
    xb <- mean(x)
    s1 <- sum((desintxb - xb)^2)
    uv <- (s1 + sum(v)) / n
    desired.sd <- sd(x)
    actualME.sd <- sqrt(uv)
    
    if (actualME.sd <= 0) stop("actualME.sd <= 0 Error")
    
    out <- desired.sd / actualME.sd
    kappa <- out - 1
    ensemble <- ensemble + kappa * (ensemble - xb)
  } else {
    kappa <- NULL
  }
  
  # Handle time series
  if (is.ts(x)) {
    ensemble <- ts(ensemble, frequency = frequency(x), start = start(x))
    dimnames(ensemble)[[2]] <- paste("Series", 1:reps)
  }
  
  list(x = x, ensemble = ensemble, xx = xx, z = z, dv = dv, dvtrim = dvtrim, 
       xmin = xmin, xmax = xmax, desintxb = desintxb, ordxx = ordxx, 
       kappa = kappa)
}