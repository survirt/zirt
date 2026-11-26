.model_label <- function(m) {
  c(grm = "Graded response model", zi = "Zero-inflated mixture graded response model",
    hurdle = "Multidimensional hurdle graded response model")[[m]]
}

#' @export
print.zirt_fit <- function(x, digits = 3, ...) {
  cat(.model_label(x$model), "\n")
  cat(sprintf("  Respondents: %d; parameters: %d; quadrature nodes: %d%s\n",
              x$n, x$npar, x$nodes,
              if (x$model == "hurdle") " per dimension" else ""))
  cat(sprintf("  logLik = %.2f; AIC = %.2f; BIC = %.2f\n", x$logLik, x$AIC, x$BIC))
  if (x$convergence != 0L) cat("  Warning: the optimiser did not report convergence.\n")
  extra <- x$coefficients[x$coefficients$parameter %in% c("pi", "rho"), ]
  if (nrow(extra) > 0L) {
    cat(sprintf("  %s = %.*f (SE %.*f, %s)\n",
                if (extra$parameter == "pi") "Proportion of structural zeros" else
                  "Correlation of susceptibility and severity",
                digits, extra$estimate, digits, extra$se, x$vcov_type))
  }
  if (!is.null(x$data$filter)) {
    cat(sprintf("  Stopped by the filter question: %.1f%% of all respondents\n",
                100 * x$data$filter$share_out))
  }
  invisible(x)
}

#' Extract item parameters
#'
#' @param object A fitted model of class \code{zirt_fit}.
#' @param ... Not used.
#' @return A data frame with one row per parameter: item, parameter name,
#'   estimate and standard error. Slopes are named \code{a} (or \code{a0}
#'   and \code{a1} for presence and severity), thresholds \code{b},
#'   \code{pi} is the proportion of structural zeros and \code{rho} the
#'   correlation of the two traits.
#' @examples
#' d <- simulate_zi(n = 800, n_items = 5, n_cat = 3, pi = 0.3, seed = 1)
#' coef(zi_fit(zirt_data(d$resp)))
#' @export
coef.zirt_fit <- function(object, ...) object$coefficients

#' @export
logLik.zirt_fit <- function(object, ...) {
  structure(object$logLik, df = object$npar, nobs = object$n, class = "logLik")
}

#' @export
summary.zirt_fit <- function(object, ...) {
  d <- object$data
  zero_obs <- sum(d$n[d$zero]) / sum(d$n)
  out <- list(model = .model_label(object$model), coefficients = object$coefficients,
              fit = c(logLik = object$logLik, npar = object$npar,
                      AIC = object$AIC, BIC = object$BIC),
              zero_observed = zero_obs, zeros = NULL, filter = d$filter)
  if (object$model == "zi") {
    pi0 <- object$coefficients$estimate[object$coefficients$parameter == "pi"]
    out$zeros <- c(observed = zero_obs, structural = pi0,
                   at_risk = max(zero_obs - pi0, 0))
  }
  class(out) <- "summary.zirt_fit"
  out
}

#' @export
print.summary.zirt_fit <- function(x, digits = 3, ...) {
  cat(x$model, "\n")
  cat(sprintf("  logLik = %.2f; parameters = %d; AIC = %.2f; BIC = %.2f\n",
              x$fit[["logLik"]], as.integer(x$fit[["npar"]]), x$fit[["AIC"]],
              x$fit[["BIC"]]))
  cat(sprintf("  All-zero respondents among those analysed: %.1f%%\n",
              100 * x$zero_observed))
  if (!is.null(x$zeros)) {
    cat(sprintf("    structural zeros: %.1f%%; zeros of respondents with the trait: %.1f%%\n",
                100 * x$zeros[["structural"]], 100 * x$zeros[["at_risk"]]))
  }
  if (!is.null(x$filter)) {
    cat(sprintf("  Stopped by the filter question: %.1f%% of all respondents\n",
                100 * x$filter$share_out))
  }
  cat("\n")
  cf <- x$coefficients
  cf$estimate <- round(cf$estimate, digits)
  cf$se <- round(cf$se, digits)
  print(cf, row.names = FALSE)
  invisible(x)
}

#' @export
print.zirt_compare <- function(x, digits = 2, ...) {
  d <- as.data.frame(x)
  d$logLik <- round(d$logLik, digits)
  d$AIC <- round(d$AIC, digits)
  d$BIC <- round(d$BIC, digits)
  d$LR <- round(d$LR, digits)
  d$p <- signif(d$p, 3)
  print(d, row.names = FALSE, na.print = "")
  cat(sprintf("Lowest BIC: %s\n", d$model[which.min(x$BIC)]))
  invisible(x)
}

#' Plot observed and expected sum-score distributions
#'
#' Draws the weighted observed distribution of the sum score as bars and the
#' distribution implied by the fitted model as points.
#'
#' @param x A fitted model of class \code{zirt_fit}.
#' @param ... Further arguments passed to \code{\link[graphics]{barplot}}.
#' @return The data frame returned by \code{\link{zi_expected}}, invisibly.
#' @examples
#' d <- simulate_zi(n = 1500, n_items = 6, n_cat = 3, pi = 0.4, seed = 2)
#' plot(zi_fit(zirt_data(d$resp)))
#' @export
plot.zirt_fit <- function(x, ...) {
  e <- zi_expected(x)
  if (anyNA(e$observed)) e$observed <- 0
  mid <- graphics::barplot(e$observed, names.arg = e$score, col = "grey80",
                           border = "grey40", xlab = "Sum score",
                           ylab = "Proportion",
                           ylim = c(0, max(e$observed, e$expected) * 1.08), ...)
  graphics::points(mid, e$expected, pch = 16)
  graphics::lines(mid, e$expected)
  graphics::legend("topright", legend = c("Observed", "Model"),
                   pch = c(15, 16), col = c("grey70", "black"), bty = "n")
  invisible(e)
}
