#' Plot multivariate time series forecast or residuals
#'
#' @param x result from \code{basicf}, \code{ridge2f} or \code{varf} (multivariate time series forecast)
#' @param selected_series name of the time series selected for plotting
#' @param type_graph "pi": basic prediction intervals;
#' "dist": a distribution of predictions; "sims": the simulations
#' @param level confidence levels for prediction intervals
#' @param ... additional parameters to be passed to \code{plot} or \code{matplot}
#'
#' @export
#'
#' @examples
#'
#' require(fpp)
#'
#' fit_obj_VAR <- ahead::varf(fpp::insurance, lags = 2,
#' h = 10, level = 95)
#'
#' fit_obj_ridge2 <- ahead::ridge2f(fpp::insurance, lags = 2,
#' h = 10, level = 95)
#'
#' par(mfrow=c(2, 2))
#' plot(fit_obj_VAR, "Quotes")
#' plot(fit_obj_VAR, "TV.advert")
#' plot(fit_obj_ridge2, "Quotes")
#' plot(fit_obj_ridge2, "TV.advert")
#'
#' obj <- ahead::ridge2f(fpp::insurance, h = 10, type_pi = "blockbootstrap",
#' block_length=5, B = 10)
#' par(mfrow=c(1, 2))
#' plot(obj, selected_series = "Quotes", type_graph = "sims",
#' main = "Predictive simulation for Quotes")
#' plot(obj, selected_series = "TV.advert", type_graph = "sims",
#' main = "Predictive simulation for TV.advert")
#'
#'
#' par(mfrow=c(1, 2))
#' plot(obj, selected_series = "Quotes", type_graph = "dist",
#' main = "Predictive simulation for Quotes")
#' plot(obj, selected_series = "TV.advert", type_graph = "dist",
#' main = "Predictive simulation for TV.advert")
#'
#'
plot.mtsforecast <- function(x, selected_series,
                             type_graph = c("pi",
                                            "dist",
                                            "sims"),
                             level = 95, ...)
{
  type_graph <- match.arg(type_graph)

  if (inherits(x, 'mtsforecast'))
  {
    if (type_graph == "pi")
    {
      y <- x$x[, selected_series]
      mean_fcast <- x$mean[, selected_series]
      upper_fcast <- x$upper[, selected_series]
      lower_fcast <- x$lower[, selected_series]

      start_y <- start(y)
      frequency_y <- frequency(y)

      y_mean <- ts(c(y, mean_fcast), start = start_y,
                   frequency = frequency_y)
      y_upper <- ts(c(y, upper_fcast), start = start_y,
                    frequency = frequency_y)
      y_lower <- ts(c(y, lower_fcast), start = start_y,
                    frequency = frequency_y)

      plot(y_mean, type='l',
           main=paste0("Forecasts for ", selected_series, " (", x$method, ")"),
           ylab="", ylim = c(min(c(y_upper, y_lower)),
                             max(c(y_upper, y_lower))), ...)
      lines(y_upper, col="gray60")
      lines(y_lower, col="gray60")
      lines(y_mean)

      return()
    }

    if (type_graph == "dist")
    {
      x <- ahead::getsimulations(x, selected_series = selected_series)$series
      start_x <- start(x)
      freq_x <- frequency(x)

      qt_upper <- apply(x, 1, function(u) quantile(u, level/100))
      qt_lower <- apply(x, 1, function(u) quantile(u, 1 - level/100))
      x_summary <- apply(x, 1, function(u) summary(u))[-3, ]
      x_ci <- rbind(apply(x, 1, function(u) t.test(u)$conf.int)[1, ],
                          rowMeans(x),
                          apply(x, 1, function(u) t.test(u)$conf.int)[2, ])
      jet.colors <- colorRampPalette( c("lightyellow", "lightgreen") )
      nbcol <- 3
      color <- jet.colors(nbcol)

      abs <- as.numeric(time(x))
      y_mean <- ts(x_summary[3, ], start = start_x, frequency = freq_x)

      bands_plot(abs, y_mean,
                 ci_upper = x_summary[1, ],
                 ci_lower = x_summary[5, ],
                 col = color[1], ylim = c(min(x_summary[1,]),
                                          max(x_summary[5,])),
                 xlab = "time",
                 ylab = selected_series,
                 ...)
      bands_add(abs, y_mean, col = color[2], ci_upper = qt_upper,
                ci_lower = qt_lower)
      bands_add(abs, y_mean, col = color[3], ci_upper = x_summary[2, ],
                ci_lower = x_summary[4, ])

      return()

    }

    if (type_graph == "sims")
    {
      B <- length(x$sims)

      actual_selected_series <- x$x[, selected_series]
      series_reps <- replicate(n = B, expr = actual_selected_series)
      series_reps <- ts(series_reps,
                        start = start(actual_selected_series),
                        frequency = frequency(actual_selected_series))
      start_x <- start(x$x)
      frequency_x <- frequency(x$x)
      start_preds <- start(x$mean)

      (preds_simulations <-
          ts(
            sapply(1:B, function (i)
              x$sims[[i]][, selected_series]),
            start = start_preds,
            frequency = frequency_x
          ))

      out <- ts(
        data = rbind(series_reps, preds_simulations),
        start = start_x,
        frequency = frequency_x
      )

      time_inputs <- as.numeric(time(x$x))
      input_series <- x$x[, selected_series]

      matplot(
        x = as.numeric(time(out)),
        y = out,
        xlab = 'time',
        ylab = selected_series,
        type = 'l',
        lwd = 2,
        ...
      )
      lines(x = time_inputs, y = input_series, lwd = 2)

      return()
    }
  } else {
    stop("Method not implemented for this type of object")
  }
}


#' Plot bivariate distributions
#'
#' @param x bivariate matrix or time series
#' @param y additional series (default is \code{NULL})
#' @export
#'
#' @examples
#'
#' obj <- ahead::ridge2f(fpp::insurance, h = 10, type_pi = "blockbootstrap",
#' block_length=5, B = 10)
#'
#' \dontrun{
#' plot(obj$residuals, type_graph = "pairs")
#' }
#'
#'
plot.pairs <- function(x, y = NULL)
{
  x <- matrix(unlist(x), ncol = 2)
  if (!is.null(y))
  {
    y <- matrix(unlist(y), ncol = 2)
    if (length(x) != length(y)) stop("'x' and 'y' must have the same dimensions")
    xvar <- c(x[,1], y[,1])
    yvar <- c(x[,2], y[,2])
    nb  <- dim(x)[1]
    zvar <- as.factor(c(rep("x", nb), rep("y", nb)))
    xy <- data.frame(xvar, yvar, zvar)
  }
  else
  {
    xvar <- x[,1]
    yvar <- x[,2]
    nb  <- dim(x)[1]
    zvar <- as.factor(rep("x", nb))
    xy <- data.frame(xvar, yvar, zvar)
  }

  #placeholder plot - prints nothing at all
  empty <- ggplot()+ggplot2::geom_point(ggplot2::aes(1,1), colour="white") +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  # plot slow af

  #scatterplot of x and y variables
  scatter <- ggplot(xy, ggplot2::aes(xvar, yvar)) +
    ggplot2::geom_point(ggplot2::aes(color=zvar)) +
    ggplot2::scale_color_manual(values = c("blue", "red")) +
    ggplot2::theme(legend.position=c(1,1),legend.justification=c(1,1))

  #marginal density of x - plot on top
  plot_top <- ggplot(xy, ggplot2::aes(xvar, fill=zvar)) +
    ggplot2::geom_density(alpha=.5) +
    ggplot2::scale_fill_manual(values = c("blue", "red")) +
    ggplot2::theme(legend.position = "none")

  #marginal density of y - plot on the right
  plot_right <- ggplot(xy, ggplot2::aes(yvar, fill=zvar)) +
    ggplot2::geom_density(alpha=.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("blue", "red")) +
    ggplot2::theme(legend.position = "none")

  #arrange the plots together, with appropriate height and width for each row and column
  grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}

# plot bands # work in progress #scratchinghead
bands_plot <- function(x, y_mean, ci_upper, ci_lower, col, y.goal = NULL, goal.col = "blue", ...)
{
  plot(x = x, y = y_mean, type = 'l', ...)
  polygon(c(x, rev(x)),
          c(ci_upper, rev(ci_lower)),
          col = col, border = FALSE)
  lines(x, y_mean, lwd = 2)
  if (!is.null(y.goal))
  {
    abline(h = y.goal, lty = 2, lwd = 2, col = goal.col)
  }
}

# add bands on a plot # work in progress #scratchinghead
bands_add <- function(x, y_mean, col, ci_upper, ci_lower)
{
  polygon(c(x, rev(x)),
          c(ci_upper, rev(ci_lower)),
          col = col, border = FALSE)
  lines(x, y_mean, lwd = 2)
}
