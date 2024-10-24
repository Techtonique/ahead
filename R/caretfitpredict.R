
#' Fit univariate time series using caret ML model (for use with \code{dynrmf})
#' 
#' @param x A matrix of predictors
#' @param y A vector of responses
#' @param method The caret method to use for fitting the model
#' @param initial_window The initial window size
#' @param horizon The forecast horizon
#' @param fixed_window Whether to use a fixed window size
#' @param tune_length Length of the tuning grid
#' @param verbose Whether to print the model summary
#' @return A model object
#' @export
fit_func <- function(x, y, method="ranger", 
initial_window = 10L,  horizon = 10L, fixed_window = FALSE, 
tune_length = 5, summary_function = NULL, 
verbose=TRUE)
{
  df <- data.frame(y=y, x)
  if (is.null(summary_function))
  {
    summary_function <- caret::defaultSummary
  }
  res <- caret::train(y ~ ., data=df,
                      method = method,
                      trControl = caret::trainControl(
                        method = "timeslice",
                        initialWindow = initial_window,                        
                        horizon = horizon,
                        fixedWindow = fixed_window,
                        skip = 0,            
                        summaryFunction = summary_function
                      ),
                      verboseIter = FALSE,
                      savePredictions = "final",
                      tuneLength = tune_length)
  if (verbose)
  {
    print(res)
  }
  return(res)
}


#' Predict univariate time series using caret ML model(for use with \code{dynrmf})
#' 
#' @param x A matrix of predictors
#' @param y A vector of responses
#' @param method The caret method to use for fitting the model
#' @param verbose Whether to print the model summary
#' @return A model object
#' @export
predict_func <- function(obj, newx)
{
  colnames(newx) <- paste0("X", seq_len(ncol(newx)))
  as.numeric(caret::predict.train(object=obj, 
                                  newdata=newx, 
                                  type = "raw"))
}

