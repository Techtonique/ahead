#' @title Ridge Regression Forecast Combination
#'
#' @description Computes forecast combination weights using Ridge Regression (OLS) regression.
#'
#' @details
#' The function integrates the Ridge Regression forecast combination implementation of the
#' \emph{ForecastCombinations} package into ForecastComb.
#'
#' The results are stored in an object of class 'ForecastComb::foreccomb_res', for which separate plot and summary functions are provided.
#'
#' @param x An object of class 'foreccomb'. Contains training set (actual values + matrix of model forecasts) and optionally a test set.
#'
#' @return Returns an object of class \code{ForecastComb::foreccomb_res} with the following components:
#' \item{Method}{Returns the best-fit forecast combination method.}
#' \item{Models}{Returns the individual input models that were used for the forecast combinations.}
#' \item{Weights}{Returns the combination weights obtained by applying the combination method to the training set.}
#' \item{Intercept}{Returns the intercept of the linear regression.}
#' \item{Fitted}{Returns the fitted values of the combination method for the training set.}
#' \item{Accuracy_Train}{Returns range of summary measures of the forecast accuracy for the training set.}
#' \item{Forecasts_Test}{Returns forecasts produced by the combination method for the test set. Only returned if input included a forecast matrix for the test set.}
#' \item{Accuracy_Test}{Returns range of summary measures of the forecast accuracy for the test set. Only returned if input included a forecast matrix and a vector of actual values for the test set.}
#' \item{Input_Data}{Returns the data forwarded to the method.}
#'
#' @examples
#' obs <- rnorm(100)
#' preds <- matrix(rnorm(1000, 1), 100, 10)
#' train_o<-obs[1:80]
#' train_p<-preds[1:80,]
#' test_o<-obs[81:100]
#' test_p<-preds[81:100,]
#'
#' data<-ForecastComb::foreccomb(train_o, train_p, test_o, test_p)
#' ahead::comb_Ridge(data)
#'
#' @seealso
#' \code{\link[ForecastCombinations]{Forecast_comb}},
#' \code{\link{foreccomb}},
#' \code{\link{plot.ForecastComb::foreccomb_res}},
#' \code{\link{summary.ForecastComb::foreccomb_res}},
#' \code{\link[forecast]{accuracy}}
#'
#' @keywords models
#'
#' @import forecast
#'
#' @export
comb_Ridge <- function(x, custom_error = NULL) {
    if (!inherits(x, "foreccomb"))
        stop("Data must be class 'foreccomb'. See ?foreccomb, to bring data in correct format.", call. = FALSE)
    observed_vector <- x$Actual_Train
    prediction_matrix <- x$Forecasts_Train
    modelnames <- x$modelnames

    lin_model <- ridge(x = prediction_matrix, y = observed_vector)
    weights <- unname(lin_model$coef[,which.min(lin_model$GCV)] )
    intercept <- lin_model$ym 
    fitted <- predict(lin_model, prediction_matrix)

    accuracy_insample <- forecast::accuracy(fitted, observed_vector)
    if (!is.null(custom_error)) {
        accuracy_insample <- cbind(accuracy_insample, custom_error(as.numeric(observed_vector), as.numeric(fitted)))
        colnames(accuracy_insample)[length(accuracy_insample)] <- "Custom Error"
    }

    if (is.null(x$Forecasts_Test) && is.null(x$Actual_Test)) {
        result <- ForecastComb::foreccomb_res(method = "Ridge Regression Regression", modelnames = modelnames, weights = weights, intercept = intercept, fitted = fitted, accuracy_insample = accuracy_insample,
                                input_data = list(Actual_Train = x$Actual_Train, Forecasts_Train = x$Forecasts_Train), 
                                predict = predict.comb_Ridge)
    }

    if (is.null(x$Forecasts_Test) == FALSE) {
        newpred_matrix <- x$Forecasts_Test
        pred <- predict(lin_model, newpred_matrix)
        if (is.null(x$Actual_Test) == TRUE) {
            result <- ForecastComb::foreccomb_res(method = "Ridge Regression Regression", modelnames = modelnames, weights = weights, intercept = intercept, fitted = fitted, accuracy_insample = accuracy_insample,
                                    pred = pred, input_data = list(Actual_Train = x$Actual_Train, Forecasts_Train = x$Forecasts_Train, Forecasts_Test = x$Forecasts_Test), 
                                    predict = predict.comb_Ridge)
        } else {
            newobs_vector <- x$Actual_Test
            accuracy_outsample <- forecast::accuracy(pred, newobs_vector)
            if (!is.null(custom_error)) {
                accuracy_outsample <- cbind(accuracy_outsample, custom_error(as.numeric(newobs_vector), as.numeric(pred)))
                colnames(accuracy_outsample)[length(accuracy_outsample)] <- "Custom Error"
            }
            result <- ForecastComb::foreccomb_res(method = "Ridge Regression Regression", modelnames = modelnames, weights = weights, intercept = intercept, fitted = fitted, accuracy_insample = accuracy_insample,
                                    pred = pred, accuracy_outsample = accuracy_outsample, input_data = list(Actual_Train = x$Actual_Train, Forecasts_Train = x$Forecasts_Train, Actual_Test = x$Actual_Test,
                                                                                                            Forecasts_Test = x$Forecasts_Test), 
                                    predict = predict.comb_Ridge)
            result$lin_model <- lin_model
        }
    }
    class(result) <- c("foreccomb_res", "comb_Ridge")
    return(result)
}

#' @export
predict.comb_Ridge <- function(object, newpreds) {
  return(predict(object$lin_model, newpreds))
}