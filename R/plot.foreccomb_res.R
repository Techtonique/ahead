#' @name plot.foreccomb_res
#' @aliases plot.foreccomb_res
#'
#' @title Plot results from forecast combination model
#' @description Produces plots for the results of a forecast combination method. Either
#' an actual vs. fitted plot (\code{which = 1}) or a barplot of the combination weights
#' (\code{which = 2}).
#'
#' @param x An object of class 'foreccomb_res'.
#' @param which Type of plot: 1 = actual vs. fitted, 2 = combination weights.
#' @param ... Other arguments passing to \code{\link[graphics]{plot.default}}.
#'
#' @return A plot for the foreccomb_res class.
#'
#' @seealso
#' \code{\link[ForecastComb]{foreccomb}},
#' \code{\link[ForecastComb]{summary.foreccomb_res}}
#'
#' @author adapted from Christoph E. Weiss and Gernot R. Roetzer (ForecastComb)
#'
#' @import ggplot2
#' @importFrom graphics barplot
#'
#' @method plot foreccomb_res
#' @export
plot.foreccomb_res <- function(x, which = 1, ...) {
  check_suggested("ggplot2")
  if (!inherits(x, "foreccomb_res"))
    stop(
      "Data must be class 'foreccomb'. See ?foreccomb, to bring data in correct format.",
      call. = FALSE
    )
  method <- x$Method
  models <- x$Models
  rolling <- FALSE
  weights <- x$Weights
  if (!is.null(dim(weights))) {
    rolling <- TRUE
    weights <- try(colMeans(weights), silent = TRUE)
    if (inherits(weights, "try-error"))
    {
      weights <- colMeans(matrix(weights, nrow = 1))
    }
  }
  fit <- x$Fitted
  forec <- x$Forecasts_Test
  observed_vector <- x$Input_Data$Actual_Train
  newobs_vector <- x$Input_Data$Actual_Test
  Index <- NULL  #Hack to satisfy CRAN check.
  
  if (which == 1) {
    if (is.null(forec) & is.null(newobs_vector)) {
      cols <- c(ACTUAL = "black",
                `COMBINED (FIT)` = "#F04546")
      
      pl <- as.data.frame(matrix(NA, ncol = 3, nrow = length(observed_vector)))
      colnames(pl) <- c("Index", "Actual", "Combined_Fit")
      pl[, 1] <- 1:nrow(pl)
      pl[, 2] <- c(observed_vector)
      pl[, 3] <- fit
      
      p <- ggplot2::ggplot(data = pl, ggplot2::aes(x = Index)) + ggplot2::geom_line(ggplot2::aes(y = Actual, colour = "ACTUAL"), size = 0.5) + ggplot2::geom_line(ggplot2::aes(y = Combined_Fit, colour = "COMBINED (FIT)"),
                                                                                                                                                                  size = 0.8) + ggplot2::scale_x_continuous(breaks = round(seq(0, max(pl$Index), by = nrow(pl) /
                                                                                                                                                                                                                                 10), 0)) + ggplot2::scale_colour_manual(name = "Series", values = cols) + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = c(0.5, 0.8)))) + ggplot2::xlab("Index") + ggplot2::ylab(paste0(method, "\n Fitted Values/Forescasts")) + ggplot2::ggtitle(paste0(
                                                                                                                                                                                                                                   method,
                                                                                                                                                                                                                                   " Forecast Combination \n Actual vs. Fitted/Test Set Forecasts"
                                                                                                                                                                                                                                 )) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")) + ggplot2::theme(legend.title = ggplot2::element_text(
          colour = "black",
          size = 12,
          face = "bold"
        ))
      p
    } else {
      cols <- c(
        ACTUAL = "black",
        `COMBINED (FIT)` = "#F04546",
        `COMBINED (FORECAST)` = "#F04546"
      )
      
      if (is.null(newobs_vector) == FALSE) {
        pl <- as.data.frame(matrix(NA, ncol = 4, nrow = (
          length(observed_vector) + length(newobs_vector)
        )))
        colnames(pl) <- c("Index",
                          "Actual",
                          "Combined_Fit",
                          "Combined_Forecast")
        pl[, 1] <- 1:nrow(pl)
        pl[, 2] <- c(observed_vector, newobs_vector)
        pl[, 3] <- c(fit, rep(NA, length(forec)))
        pl[, 4] <- c(rep(NA, length(fit)), forec)
        pl[length(observed_vector), 4] <- pl[length(observed_vector), 3]
        
        p <- ggplot2::ggplot(data = pl, ggplot2::aes(x = Index)) + ggplot2::geom_line(ggplot2::aes(y = Actual, colour = "ACTUAL"), size = 0.5) + ggplot2::geom_line(ggplot2::aes(y = c(Combined_Fit), colour = "COMBINED (FIT)"),
                                                                                                                                                                    na.rm = TRUE,
                                                                                                                                                                    size = 0.8) + ggplot2::geom_line(ggplot2::aes(y = c(Combined_Forecast), colour = "COMBINED (FORECAST)"),
                                                                                                                                                                                                     na.rm = TRUE,
                                                                                                                                                                                                     size = 1.5) + ggplot2::scale_x_continuous(breaks = round(seq(0, max(pl$Index), by = nrow(pl) /
                                                                                                                                                                                                                                                                    10), 0)) + ggplot2::scale_colour_manual(name = "Series", values = cols) + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = c(
                                                                                                                                                                                                                                                                      0.5, 0.8, 1.5
                                                                                                                                                                                                                                                                    )))) + ggplot2::xlab("Index") + ggplot2::ylab(paste0(method, "\n Fitted Values/Forecasts")) + ggplot2::ggtitle(paste0(
                                                                                                                                                                                                                                                                      method,
                                                                                                                                                                                                                                                                      " Forecast Combination \n Actual vs. Fitted/Test Set Forecast"
                                                                                                                                                                                                                                                                    )) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")) + ggplot2::theme(legend.title = ggplot2::element_text(
            colour = "black",
            size = 12,
            face = "bold"
          )) +
          ggplot2::geom_vline(
            xintercept = length(observed_vector),
            size = 1,
            linetype = "longdash",
            colour = "black"
          )
        p
      } else {
        pl <- as.data.frame(matrix(NA, ncol = 4, nrow = (
          length(observed_vector) + length(forec)
        )))
        colnames(pl) <- c("Index",
                          "Actual",
                          "Combined_Fit",
                          "Combined_Forecast")
        pl[, 1] <- 1:nrow(pl)
        pl[, 2] <- c(observed_vector, rep(NA, length(forec)))
        pl[, 3] <- c(fit, rep(NA, length(forec)))
        pl[, 4] <- c(rep(NA, length(fit)), forec)
        pl[length(observed_vector), 4] <- pl[length(observed_vector), 3]
        
        p <- ggplot2::ggplot(data = pl, ggplot2::aes(x = Index)) + ggplot2::geom_line(ggplot2::aes(y = Actual, colour = "ACTUAL"),
                                                                                      na.rm = TRUE,
                                                                                      size = 0.5) + ggplot2::geom_line(ggplot2::aes(y = c(Combined_Fit), colour = "COMBINED (FIT)"),
                                                                                                                       na.rm = TRUE,
                                                                                                                       size = 0.8) + ggplot2::geom_line(ggplot2::aes(y = c(Combined_Forecast), colour = "COMBINED (FORECAST)"),
                                                                                                                                                        na.rm = TRUE,
                                                                                                                                                        size = 1.5) +
          ggplot2::scale_x_continuous(breaks = round(seq(0, max(pl$Index), by = nrow(pl) /
                                                           10), 0)) + ggplot2::scale_colour_manual(name = "Series", values = cols) + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = c(
                                                             0.5, 0.8, 1.5
                                                           )))) + ggplot2::xlab("Index") + ggplot2::ylab(paste0(method, "\n Fitted Values")) + ggplot2::ggtitle(paste0(method, " Forecast Combination \n Actual vs. Fitted")) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")) + ggplot2::theme(legend.title = ggplot2::element_text(
            colour = "black",
            size = 12,
            face = "bold"
          )) + ggplot2::geom_vline(
            xintercept = length(observed_vector),
            size = 1,
            linetype = "longdash",
            colour = "black"
          )
        p
      }
    }
  } else {
    if (which == 2) {
      if (is.numeric(weights)) {
        if (!rolling) {
          graphics::barplot(
            weights,
            main = paste0(method, "\nCombination Weights"),
            ylab = "Combination Weight",
            names.arg = models,
            ylim = c(min(1.1 * min(weights), 0), 1.1 * max(weights)),
            las = 3,
            cex.names = 0.8
          )
        } else{
          graphics::barplot(
            weights,
            main = paste0(method, "\n Average Rolling Combination Weights"),
            ylab = "Average Combination Weight",
            names.arg = models,
            ylim = c(min(1.1 * min(weights), 0), 1.1 * max(weights)),
            las = 3,
            cex.names = 0.8
          )
        }
      } else {
        message(
          paste0(
            method,
            " produces time-varying weights among input models. Cannot plot weights."
          )
        )
      }
    } else
      stop("Parameter 'which' must be either 1 or 2.", call. = FALSE)
  }
}