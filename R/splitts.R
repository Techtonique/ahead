#' Partition a time series object
#' 
#' @param y A time series object
#' @param split_prob Splitting ratio 
#' @param return_indices if TRUE, returns series' indices, otherwise, time series objects
#' 
#' @export
#' @examples
#' 
#' misc::splitts(ts(1:10))
#' 
splitts <-
  function(y,
           split_prob = 0.5,
           return_indices = FALSE)
  {
    n_y <- base::ifelse(test = is.null(dim(y)),
                        yes = length(y),
                        no = dim(y)[1])

    index_train <- 1:floor(split_prob * n_y)
    if (return_indices)
      return(index_train)

    start_y <- stats::start(y)
    frequency_y <- stats::frequency(y)

    if (is.null(ncol(y)))
      # univariate case
    {
      training <- ts(y[index_train],
                     start = start_y,
                     frequency = frequency_y)
      start_testing <- tsp(training)[2] + 1 / frequency_y
      return(list(
        training = training,
        testing = ts(y[-index_train],
                     start = start_testing,
                     frequency = frequency_y)
      ))
    } else {
      # multivariate case
      training <- ts(y[index_train, ],
                     start = start_y,
                     frequency = frequency_y)
      start_testing <- tsp(training)[2] + 1 / frequency_y
      return(list(
        training = training,
        testing = ts(y[-index_train, ],
                     start = start_testing,
                     frequency = frequency_y)
      ))
    }
  }