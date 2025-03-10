#' Compute global attention weights and context vectors for time series
#' 
#' @param series Numeric vector containing the time series of length n
#' @return List containing:
#'   \item{attention_weights}{n × n matrix where entry (i,j) represents the attention 
#'         weight of time j on time i. Only entries j <= i are non-zero (causal attention).}
#'   \item{context_vectors}{Vector of length n where each entry i is the weighted sum 
#'         of all values up to time i, using the attention weights.}
#' @examples
#' # For a series of length 5
#' series <- c(1, 2, 3, 4, 5)
#' result <- compute_attention(series)
#' 
#' # attention_weights will be 5x5 matrix
#' # context_vectors will be length 5
#' dim(result$attention_weights)  # [1] 5 5
#' length(result$context_vectors) # [1] 5
#' @export
computeattention <- function(series) {
  # Input validation
  if (!is.numeric(series)) stop("series must be numeric")
  
  # Compute attention weights using Rcpp function
  attention_weights <- compute_attention_cpp(series)
  
  # Compute context vectors using Rcpp function
  context_vectors <- compute_context_vectors_cpp(series, attention_weights)
  
  # Add dimension names for clarity
  n <- length(series)
  dimnames(attention_weights) <- list(
    paste0("t", 1:n),
    paste0("t", 1:n)
  )
  names(context_vectors) <- paste0("t", 1:n)
  
  return(list(
    attention_weights = attention_weights,
    context_vectors = context_vectors
  ))
} 