#' Compute global attention weights and context vectors for time series
#'
#' @param series Numeric vector containing the time series of length n
#' @param attention_type String specifying the type of attention mechanism to use.
#'        Options are: "cosine", "exponential", "dot_product", "scaled_dot_product",
#'        "gaussian", "linear", "value_based", "hybrid", "parametric". Default is "cosine".
#' @param window_size Integer parameter for window size (applicable for "cosine" attention).
#' @param decay_factor Double for decay factor (applicable for "exponential" attention).
#' @param temperature Double for temperature (applicable for "scaled_dot_product" attention).
#' @param sigma Double for sigma (applicable for "gaussian" attention).
#' @param sensitivity Double for sensitivity (applicable for "value_based" or "hybrid" attention).
#' @param alpha Double for alpha (applicable for "parametric" attention).
#' @param beta Double for beta (applicable for "parametric" attention).
#'
#' @return List containing:
#'   \item{attention_weights}{n × n matrix where entry (i,j) represents the attention
#'         weight of time j on time i. Only entries j <= i are non-zero (causal attention).}
#'   \item{context_vectors}{Vector of length n where each entry i is the weighted sum
#'         of all values up to time i, using the attention weights.}
#'
#' @examples
#' # For a series of length 5 using "cosine" attention
#' series <- c(1, 2, 3, 4, 5)
#' result <- computeattention(series, attention_type = "cosine", window_size = 3)
#'
#' # attention_weights will be 5x5 matrix
#' # context_vectors will be length 5
#' dim(result$attention_weights)  # [1] 5 5
#' length(result$context_vectors) # [1] 5
#'
#' @export
computeattention <- function(series,
                             attention_type = "cosine",
                             window_size = 3,
                             decay_factor = 5.0,
                             temperature = 1.0,
                             sigma = 1.0,
                             sensitivity = 1.0,
                             alpha = 0.5,
                             beta = 0.5) {

  # Input validation
  if (!is.numeric(series)) stop("series must be numeric")

  # Check if the attention_type is valid
  valid_attention_types <- c("cosine", "exponential", "dot_product",
                             "scaled_dot_product", "gaussian", "linear",
                             "value_based", "hybrid", "parametric")
  if (!(attention_type %in% valid_attention_types)) {
    stop(paste("Invalid attention_type. Choose from:", paste(valid_attention_types, collapse = ", ")))
  }

  # Compute attention weights based on the chosen attention mechanism
  attention_weights <- switch(attention_type,
                              "cosine" = cosine_attention_cpp(series, window_size),
                              "exponential" = exponential_attention_cpp(series, decay_factor),
                              "dot_product" = dot_product_attention_cpp(series),
                              "scaled_dot_product" = scaled_dot_product_attention_cpp(series, temperature),
                              "gaussian" = gaussian_attention_cpp(series, sigma),
                              "linear" = linear_attention_cpp(series),
                              "value_based" = value_based_attention_cpp(series, sensitivity),
                              "hybrid" = hybrid_attention_cpp(series, decay_factor, sensitivity),
                              "parametric" = parametric_attention_cpp(series, alpha, beta)
  )

  # Compute context vectors using the attention weights
  context_vectors <- compute_context_vectors_cpp(series, attention_weights)

  # Add dimension names for clarity
  n <- length(series)
  dimnames(attention_weights) <- list(
    paste0("t", 1:n),
    paste0("t", 1:n)
  )
  names(context_vectors) <- paste0("t", 1:n)

  # Return a list with attention weights and context vectors
  return(list(
    attention_weights = attention_weights,
    context_vectors = context_vectors
  ))
}
