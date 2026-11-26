#' Generate responses with structural zeros or a hurdle process
#'
#' Simulates ordinal item responses from the zero-inflated mixture graded
#' response model or from the multidimensional hurdle graded response model,
#' optionally preceded by a scale-level filter question.
#'
#' @details
#' With \code{model = "zi"} a proportion \code{pi} of respondents are
#' structural zeros and the others answer according to a graded response
#' model. With \code{model = "hurdle"} each item is left at zero unless a
#' presence process in the susceptibility trait is passed, and the nonzero
#' category depends on the severity trait; the two traits correlate
#' \code{rho}. Default thresholds place the scale at the upper end of the
#' trait, as is typical of symptom and problem scales.
#'
#' If \code{filter_sens} is given, a filter question is generated: every
#' respondent with a nonzero pattern passes it, respondents with an all-zero
#' pattern from the trait class pass with probability \code{filter_sens},
#' and structural zeros do not pass. Items of respondents who do not pass
#' are set to \code{NA}.
#'
#' @param n Number of respondents.
#' @param n_items Number of items.
#' @param n_cat Number of response categories per item (2 for binary).
#' @param pi Proportion of structural zeros (zero-inflated model).
#' @param model \code{"zi"} or \code{"hurdle"}.
#' @param rho Correlation of susceptibility and severity (hurdle model).
#' @param a Optional vector of slopes.
#' @param b Optional list with one vector of thresholds per item (for the
#'   hurdle model, the susceptibility threshold first and the severity
#'   thresholds after it).
#' @param filter_sens Optional probability that an at-risk respondent with
#'   an all-zero pattern passes the filter question.
#' @param seed Optional random seed.
#'
#' @return A list with \code{resp}, \code{structural} (indicator of the
#'   structural-zero class), the traits, the generating parameters and, if
#'   requested, \code{filter}.
#'
#' @examples
#' d <- simulate_zi(n = 500, n_items = 5, n_cat = 4, pi = 0.3, seed = 1)
#' table(rowSums(d$resp) == 0, d$structural)
#' @export
simulate_zi <- function(n = 1000L, n_items = 6L, n_cat = 4L, pi = 0.3,
                        model = c("zi", "hurdle"), rho = 0.5, a = NULL,
                        b = NULL, filter_sens = NULL, seed = NULL) {
  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)
  if (n_cat < 2L) stop("'n_cat' must be at least 2.", call. = FALSE)
  if (is.null(a)) a <- seq(1.2, 2.2, length.out = n_items)
  if (is.null(b)) {
    first <- seq(0, 1, length.out = n_items)
    b <- lapply(first, function(f) f + 0.8 * (seq_len(n_cat - 1L) - 1L))
    if (model == "hurdle") {
      b <- lapply(first, function(f) c(f, -1 + 1.0 * (seq_len(max(n_cat - 2L, 0L)) - 1L)))
    }
  }
  draw <- function(P) {
    u <- stats::runif(nrow(P))
    rowSums(u > t(apply(P, 1L, cumsum))[, -ncol(P), drop = FALSE])
  }
  resp <- matrix(0L, n, n_items)
  if (model == "zi") {
    structural <- stats::rbinom(n, 1L, pi)
    theta <- stats::rnorm(n)
    for (j in seq_len(n_items)) {
      resp[, j] <- draw(p_grm(theta, a[j], b[[j]]))
    }
    resp[structural == 1L, ] <- 0L
    traits <- list(theta = theta)
  } else {
    structural <- rep(0L, n)
    t0 <- stats::rnorm(n)
    t1 <- rho * t0 + sqrt(1 - rho^2) * stats::rnorm(n)
    for (j in seq_len(n_items)) {
      resp[, j] <- draw(p_hurdle(t0, t1, a[j], b[[j]][1L], a[j], b[[j]][-1L]))
    }
    traits <- list(susceptibility = t0, severity = t1)
  }
  storage.mode(resp) <- "integer"
  colnames(resp) <- sprintf("y%d", seq_len(n_items))
  out <- c(list(resp = resp, structural = structural), traits,
           list(a = a, b = b, pi = if (model == "zi") pi else NA_real_,
                rho = if (model == "hurdle") rho else NA_real_))
  if (!is.null(filter_sens)) {
    zero <- rowSums(resp) == 0L
    pass <- ifelse(!zero, 1L, ifelse(structural == 1L, 0L,
                                     stats::rbinom(n, 1L, filter_sens)))
    resp[pass == 0L, ] <- NA_integer_
    out$resp <- resp
    out$filter <- pass
  }
  out
}
