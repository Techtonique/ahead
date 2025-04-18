#' @name summary.foreccomb_res
#' @aliases print.foreccomb_res_summary

#' @title Summary of Forecast Combination
#' @description \code{summary} method for class \sQuote{foreccomb_res}. Includes information about combination method,
#' combination weights assigned to the individual forecast models, as well as an accuracy evaluation of the combined
#' forecast.
#'
#' @param object An object of class 'foreccomb'. Contains training set (actual values + matrix of model forecasts) and optionally a test set.
#' @param x An object of class 'foreccomb'. Contains training set (actual values + matrix of model forecasts) and optionally a test set.
#' @param ... potential further arguments (require by generic)
#'
#'
#' @seealso
#' \code{\link[ForecastComb]{foreccomb}},
#' \code{\link[ForecastComb]{plot.foreccomb_res}},
#'
#' @author Christoph E. Weiss and Gernot R. Roetzer
#'
#' @import forecast
#'
#' @rdname summary.foreccomb_res
#' @method summary foreccomb_res
#' @export
summary.foreccomb_res <- function(object, ...) {
    if (!inherits(object, "foreccomb_res"))
        stop("Data must be class 'foreccomb'. See ?foreccomb, to bring data in correct format.", call. = FALSE)

    ans <- list()

    ans$Method <- object$Method

    if (!is.character(object$Weights)) {
      if(is.null(dim(object$Weights))){
        ans$weight <- matrix(object$Weights, ncol = 1)
        colnames(ans$weight) <- "Combination Weight"
        rownames(ans$weight) <- object$Models
      }
      else{
        ans$weight <- matrix(object$Weights[dim(object$Weights)[1],], ncol = 1)
        colnames(ans$weight) <- "Combination Weight (End of Test Set)"
        rownames(ans$weight) <- object$Models
      }
    } else {
        ans$weight <- "Weights of the individual forecasts differ over time with trimmed mean"
    }

    ans$Intercept <- object$Intercept
    if(!is.null(ans$Intercept) & length(ans$Intercept)>1)
      ans$Intercept <- ans$Intercept[length(ans$Intercept)]

    ans$accuracy <- as.data.frame(rbind(object$Accuracy_Train[1:5], object$Accuracy_Test))
    rownames(ans$accuracy)[1] <- "Training Set"
    colnames(ans$accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")

    ans$data <- deparse(substitute(object))

    ans <- append(ans, subset(object, !(names(object) %in% c("Method", "Weights", "Intercept", "Accuracy_Train", "Accuracy_Test"))))

    class(ans) <- c("foreccomb_res_summary")

    return(ans)
}

#' @rdname summary.foreccomb_res
#' @method print foreccomb_res_summary
#' @export
print.foreccomb_res_summary <- function(x, ...) {
    if (inherits(x, "foreccomb_res_summary"))
        stop("Data must be class 'foreccomb_res_summary'", call. = FALSE)

    cat("\n")
    cat("Summary of Forecast Combination \n")
    cat("------------------------------- \n")
    cat("\n")
    cat("Method: ", x$Method, "\n")
    cat("\n")
    cat("Individual Forecasts & Combination Weights: \n")
    if (!is.character(x$weight)) {
        cat("\n")
        print(x$weight)
    } else {
        cat(x$weight)
        cat("\n")
    }
    cat("\n")
    if (!is.null(x$Intercept)) {
        cat("Intercept (Bias-Correction): ", x$Intercept, "\n")
        cat("\n")
    }
    cat("Accuracy of Combined Forecast: \n")
    cat("\n")
    print(x$accuracy)
    cat("\n")
    cat("Additional information can be extracted from the combination object: \n")
    cat("For fitted values (training set): ", paste0(x$data, "$Fitted"), "\n")
    if (!is.null(x$Forecasts_Test)) {
        cat("For forecasts (test set): ", paste0(x$data, "$Forecasts_Test"), "\n")
    }
    cat("See ", paste0("str(", x$data, ")"), " for full list.")
}