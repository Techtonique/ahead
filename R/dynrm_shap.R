# dynrmf_shap.R
# Exact Shapley values for ahead::dynrmf.
# Mirrors the dynrmf() argument list — no manual model-list needed.

#' Compute exact Shapley values for an ahead::dynrmf model
#'
#' @param y            Training time series (passed to ahead::dynrmf as y)
#' @param xreg_fit     Regressor matrix used during training
#' @param xreg_predict Regressor matrix for the forecast horizon (h x p);
#'                     h and feature names are derived from this matrix.
#' @param fit_func     Fitting function for ahead::dynrmf (default: ahead::ridge)
#' @param verbose      Print progress per horizon and efficiency check (default: TRUE)
#' @param tol          Tolerance for the efficiency-axiom check (default: 1e-10)
#' @export
#' @return A list of class "dynrmf_shap" with:
#'   $S           : h x p matrix of Shapley values (per variable, per horizon)
#'   $X           : test regressor matrix
#'   $baseline    : scalar baseline (mean of empty-coalition forecasts)
#'   $baselines   : per-horizon baseline forecasts
#'   $predictions : per-horizon full-coalition forecasts
#'   $feature_names
#'   $h
dynrmf_shap <- function(y,
                            xreg_fit,
                            xreg_predict,
                            fit_func = ahead::ridge,
                            verbose  = TRUE,
                            tol      = 1e-10) {
  
  # ── Coerce inputs ──────────────────────────────────────────────────────────
  if (is.ts(xreg_predict) || is.data.frame(xreg_predict))
    xreg_predict <- as.matrix(xreg_predict)
  if (is.ts(xreg_fit) || is.data.frame(xreg_fit))
    xreg_fit <- as.matrix(xreg_fit)
  
  # ── Derive h, p, feature_names from xreg_predict ──────────────────────────
  h             <- nrow(xreg_predict)
  p             <- ncol(xreg_predict)
  feature_names <- if (!is.null(colnames(xreg_predict))) colnames(xreg_predict) else
    paste0("x", seq_len(p))
  colnames(xreg_predict) <- feature_names
  
  # ── Reference values: training-set column means ───────────────────────────
  ref        <- colMeans(xreg_fit)
  names(ref) <- feature_names
  
  # ── Helper: call dynrmf with a constant regressor matrix, return mean[tau] ─
  .v <- function(coal_vec, tau) {
    xmat <- matrix(rep(coal_vec, h), nrow = h, ncol = p, byrow = TRUE)
    colnames(xmat) <- feature_names
    ahead::dynrmf(
      y            = y,
      xreg_fit     = xreg_fit,
      xreg_predict = xmat,
      fit_func     = fit_func,
      h            = h
    )$mean[tau]
  }
  
  # ── Coalition enumeration (computed once, outside the tau loop) ────────────
  n_coal <- 2L ^ p
  masks  <- 0L:(n_coal - 1L)
  jbits  <- bitwShiftL(1L, seq_len(p) - 1L)   # jbits[j] = bitmask for feature j
  
  # membership[k, j] = TRUE  <=>  feature j is in coalition masks[k]
  membership <- outer(masks, jbits, function(m, b) bitwAnd(m, b) > 0L)
  
  # Coalition sizes 0 … p
  s_sizes <- rowSums(membership)
  
  # Shapley weights by coalition size (log-factorials avoid overflow for large p)
  log_fp    <- lfactorial(p)
  w_by_size <- exp(
    lfactorial(seq_len(p) - 1L) +   # lfactorial(s),     s = 0, …, p-1
      lfactorial(p - seq_len(p))   -   # lfactorial(p-s-1)
      log_fp
  )
  
  # ── Main loop over forecast horizons ──────────────────────────────────────
  Phi            <- matrix(0.0, nrow = h, ncol = p,
                           dimnames = list(NULL, feature_names))
  baseline_vec   <- numeric(h)
  prediction_vec <- numeric(h)
  
  for (tau in seq_len(h)) {
    if (verbose) message(sprintf("[dynrmf_shap] horizon %d / %d", tau, h))
    
    # Cache v(S) for every coalition at this horizon
    v <- numeric(n_coal)
    for (k in seq_len(n_coal)) {
      coal <- ref
      active_cols <- which(membership[k, ])
      if (length(active_cols) > 0L)
        coal[active_cols] <- xreg_predict[tau, active_cols]
      v[k] <- .v(coal, tau)
    }
    
    baseline_vec[tau]   <- v[1L]       # empty coalition  (mask 0)
    prediction_vec[tau] <- v[n_coal]   # full coalition   (mask 2^p - 1)
    
    # Vectorised Shapley summation — no inner loop over coalitions
    for (j in seq_len(p)) {
      not_j    <- !membership[, j]            # 2^(p-1) coalitions without j
      m_no_j   <- masks[not_j]
      m_with_j <- bitwOr(m_no_j, jbits[j])
      
      weights     <- w_by_size[s_sizes[not_j] + 1L]
      Phi[tau, j] <- sum(weights * (v[m_with_j + 1L] - v[m_no_j + 1L]))
    }
  }
  
  # ── Checks ────────────────────────────────────────────────────────────────
  # 1. Efficiency axiom: rowSums(Phi) == predictions - baselines, per horizon
  residuals <- rowSums(Phi) - (prediction_vec - baseline_vec)
  max_resid <- max(abs(residuals))
  axiom_ok  <- max_resid <= tol
  
  # 2. Total Shapley sum == sum(predictions) - sum(baselines)
  total_gap     <- sum(prediction_vec) - sum(baseline_vec)
  total_diff    <- abs(sum(Phi) - total_gap)
  total_shap_ok <- total_diff <= tol
  
  # 3. scalar baseline is mean(baselines) — always true by construction, shown for transparency
  baseline_scalar <- mean(baseline_vec)
  
  if (verbose) {
    cat(sprintf(
      "[dynrmf_shap] (1) Efficiency axiom   max|rowSums(S) - (pred - base)| = %.3e  [%s]\n",
      max_resid, if (axiom_ok) "PASS" else "FAIL"
    ))
    cat(sprintf(
      "[dynrmf_shap] (2) Total Shapley sum  |sum(S) - (sum(pred) - sum(base))| = %.3e  [%s]\n",
      total_diff, if (total_shap_ok) "PASS" else "FAIL"
    ))
    cat(sprintf(
      "[dynrmf_shap] (3) baseline == mean(baselines): %.6f == %.6f  [%s]\n",
      baseline_scalar, mean(baseline_vec), "PASS"
    ))
  }
  
  if (!axiom_ok)
    warning(sprintf(
      "Efficiency axiom violated: max|rowSums(S) - (predictions - baselines)| = %.3e (tol = %.3e).",
      max_resid, tol
    ))
  if (!total_shap_ok)
    warning(sprintf(
      "Total Shapley sum check failed: |sum(S) - (sum(pred) - sum(base))| = %.3e (tol = %.3e).",
      total_diff, tol
    ))
  
  structure(
    list(
      S             = Phi,
      X             = xreg_predict,
      baseline      = baseline_scalar,
      baselines     = baseline_vec,
      predictions   = prediction_vec,
      feature_names = feature_names,
      h             = h
    ),
    class = "dynrmf_shap"
  )
}

#' Waterfall plot for a dynrmf_shap object
#'
#' Aggregates SHAP values across horizons and draws a horizontal waterfall.
#'
#' @param shap     A "dynrmf_shap" object from dynrmf_shap()
#' @param title    Plot title string
#' @param agg_fun  How to aggregate across horizons: sum (default) or mean
#' @param colors   Named list with elements "pos", "neg", "base"
#' @export
#' @return A ggplot2 object
plot_dynrmf_shap_waterfall <- function(shap,
                                title   = "SHAP Waterfall",
                                agg_fun = sum,
                                colors  = list(pos  = "#4575b4",
                                               neg  = "#d73027",
                                               base = "#878787")) {
  
  stopifnot(inherits(shap, "dynrmf_shap"))
  
  phi  <- apply(shap$S, 2, agg_fun)       # aggregate Shapley per feature
  f_x  <- agg_fun(shap$predictions)       # aggregated model output
  E_x  <- agg_fun(shap$baselines)         # aggregated baseline
  
  # Sort by |phi| ascending so largest bar is at the top
  ord  <- order(abs(phi))
  phi  <- phi[ord]
  nms  <- names(phi)
  n    <- length(phi)
  
  # Running x positions for waterfall bars
  starts <- E_x + c(0, cumsum(phi[-n]))
  ends   <- E_x + cumsum(phi)
  
  feat_levels <- c("E[f(X)]", nms, "f(x)")
  
  # Feature bars
  df <- data.frame(
    feature  = factor(nms, levels = feat_levels),
    xbeg     = starts,
    xend     = ends,
    phi      = phi,
    bar_type = ifelse(phi >= 0, "pos", "neg"),
    stringsAsFactors = FALSE
  )
  
  # Baseline row
  base_row <- data.frame(
    feature  = factor("E[f(X)]", levels = feat_levels),
    xbeg     = 0,   xend = E_x,
    phi      = E_x, bar_type = "base",
    stringsAsFactors = FALSE
  )
  
  # Total row
  total_row <- data.frame(
    feature  = factor("f(x)", levels = feat_levels),
    xbeg     = 0,   xend = f_x,
    phi      = f_x, bar_type = "base",
    stringsAsFactors = FALSE
  )
  
  all_df <- rbind(base_row, df, total_row)
  
  # Label: signed contribution for feature bars, raw value for end-points
  x_range  <- diff(range(c(all_df$xbeg, all_df$xend)))
  pad      <- 0.012 * x_range
  all_df$label   <- sprintf("%+.3f", all_df$phi)
  is_endpoint     <- all_df$feature %in% c("E[f(X)]", "f(x)")
  all_df$label[is_endpoint] <- sprintf("%.3f", all_df$xend[is_endpoint])
  all_df$label_x  <- pmax(all_df$xbeg, all_df$xend) + pad
  
  col_map <- c(pos = colors$pos, neg = colors$neg, base = colors$base)
  
  ggplot2::ggplot(all_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        ymin  = as.numeric(feature) - 0.38,
        ymax  = as.numeric(feature) + 0.38,
        xmin  = xbeg,
        xmax  = xend,
        fill  = bar_type
      ),
      colour    = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(
      xintercept = E_x,
      linetype   = "dashed",
      colour     = "grey45",
      linewidth  = 0.4
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = label_x, y = feature, label = label),
      hjust = 0,
      size  = 3.0
    ) +
    ggplot2::scale_fill_manual(values = col_map, guide = "none") +
    ggplot2::scale_y_discrete(limits = feat_levels) +
    ggplot2::expand_limits(x = max(all_df$label_x) + pad * 3) +
    ggplot2::labs(title = title, x = "Aggregated contribution", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold", size = 11)
    )
}



# =============================================================================
# USAGE EXAMPLE
# =============================================================================
#
# library(fpp2); library(ahead); library(e1071); library(misc)
# library(ggplot2); library(patchwork)
# 
# y       <- fpp2::uschange[, "Consumption"]
# xreg    <- scale(fpp2::uschange[, c("Income", "Savings", "Unemployment")])
# split   <- misc::splitts(y, split_prob = 0.9)
# xreg_train <- window(xreg, start = start(split$training), end = end(split$training))
# xreg_test <- window(xreg, start = start(split$testing),  end = end(split$testing))
# 
# shap <- ahead::dynrmf_shap(
#   y            = split$training,
#   xreg_fit     = xreg_train,
#   xreg_predict = xreg_test,
#   fit_func     = e1071::svm
# )
# 
# p1 <- ahead::plot_dynrmf_shap_waterfall(shap, title = "Baseline scenario")
# 
# xreg_pess <- xreg_test
# xreg_pess[,"Income"] <- -1;  
# xreg_pess[,"Savings"] <- -0.5
# 
# shap_pess <- dynrmf_shap(
#   y            = split$training,
#   xreg_fit     = xreg_train,
#   xreg_predict = xreg_pess,
#   fit_func     = e1071::svm
# )
# 
# p2 <- ahead::plot_dynrmf_shap_waterfall(shap_pess, title = "Pessimistic scenario")
# 
# xreg_opt  <- xreg_test
# xreg_opt[,"Income"]  <-  2;  
# xreg_opt[,"Savings"]  <-  0.5
# 
# shap_opt <- dynrmf_shap(
#   y            = split$training,
#   xreg_fit     = xreg_train,
#   xreg_predict = xreg_opt,
#   fit_func     = e1071::svm
# )
# 
# p3 <- plot_dynrmf_shap_waterfall(shap_opt, title = "Optimistic scenario")
# 
# xreg_ovr  <- xreg_test
# xreg_ovr[,"Income"]  <-  2.5; 
# xreg_ovr[,"Savings"] <-  0.75
# 
# shap_ovr <- dynrmf_shap(
#   y            = split$training,
#   xreg_fit     = xreg_train,
#   xreg_predict = xreg_ovr,
#   fit_func     = e1071::svm
# )
# 
# p4 <- plot_dynrmf_shap_waterfall(shap_ovr, title = "Overly optimistic scenario")
# 
# # (p1 + p2)/(p3 + p4)
