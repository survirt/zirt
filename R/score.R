# Posterior weights of the nodes for every pattern, and node positions.
.posterior <- function(fit) {
  d <- fit$data
  p <- split_par(fit$par, .sizes(fit$model, d$K))
  nodes <- if (fit$model == "hurdle") gh_grid2(fit$nodes, tanh(p$zrho)) else
    gh_nodes(fit$nodes)
  L <- pattern_lik(d$pat, .plist(p, fit$model, d$K, nodes))
  Lw <- sweep(L, 2L, nodes$w, `*`)
  list(Lw = Lw, nodes = nodes, p = p)
}

#' Score respondents under a fitted model
#'
#' Computes expected a posteriori (EAP) trait estimates with posterior
#' standard deviations. For the zero-inflated model it also returns the
#' posterior probability that a respondent is a structural zero; for the
#' hurdle model it returns susceptibility and severity scores.
#'
#' @details
#' Under the zero-inflated model a respondent with an all-zero pattern is a
#' structural zero with probability
#' \eqn{\pi / \{\pi + (1 - \pi) P(y = 0 \mid \mathrm{trait\ class})\}}, and
#' every other respondent has probability 0. The trait score of the
#' zero-inflated model refers to the class that has the trait.
#'
#' @param fit A fitted model of class \code{zirt_fit}.
#' @param patterns If \code{TRUE}, return one row per distinct response
#'   pattern with its weighted frequency instead of one row per respondent.
#'
#' @return A data frame. Columns are \code{theta} and \code{sd} (graded
#'   response and zero-inflated models), \code{p_structural} (zero-inflated
#'   model), or \code{susceptibility}, \code{sd_susceptibility},
#'   \code{severity} and \code{sd_severity} (hurdle model).
#'
#' @examples
#' d <- simulate_zi(n = 1500, n_items = 6, n_cat = 3, pi = 0.4, seed = 2)
#' fit <- zi_fit(zirt_data(d$resp))
#' head(zi_score(fit, patterns = TRUE))
#' @export
zi_score <- function(fit, patterns = FALSE) {
  stopifnot(inherits(fit, "zirt_fit"))
  d <- fit$data
  po <- .posterior(fit)
  marg <- rowSums(po$Lw)
  W <- po$Lw / marg
  eap <- function(x) {
    m <- drop(W %*% x)
    cbind(m, sqrt(pmax(drop(W %*% x^2) - m^2, 0)))
  }
  if (fit$model == "hurdle") {
    s0 <- eap(po$nodes$t0)
    s1 <- eap(po$nodes$t1)
    out <- data.frame(susceptibility = s0[, 1L], sd_susceptibility = s0[, 2L],
                      severity = s1[, 1L], sd_severity = s1[, 2L])
  } else {
    s <- eap(po$nodes$x)
    out <- data.frame(theta = s[, 1L], sd = s[, 2L])
    if (fit$model == "zi") {
      pi0 <- stats::plogis(po$p$logitpi)
      out$p_structural <- ifelse(d$zero, pi0 / (pi0 + (1 - pi0) * marg), 0)
    }
  }
  if (patterns) {
    out <- cbind(as.data.frame(d$pat), n = d$n, out)
    names(out)[seq_along(d$items)] <- d$items
    return(out)
  }
  out <- out[d$idx, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Distribution of the sum score at each node (generalised Lord-Wingersky).
.sum_dist <- function(plist) {
  nn <- nrow(plist[[1L]])
  dist <- matrix(1, nn, 1L)
  for (P in plist) {
    K <- ncol(P)
    new <- matrix(0, nn, ncol(dist) + K - 1L)
    for (k in seq_len(K)) {
      cols <- seq_len(ncol(dist)) + k - 1L
      new[, cols] <- new[, cols] + dist * P[, k]
    }
    dist <- new
  }
  dist
}

#' Observed and model-implied sum-score distributions
#'
#' Returns the weighted observed distribution of the sum score next to the
#' distribution implied by a fitted model. The comparison shows how well a
#' model reproduces the pile of zero scores and the shape of the nonzero
#' scores. Only respondents with complete responses enter the observed
#' distribution. When more than 5\% of the weighted sample has missing
#' responses, as happens when follow-up items are asked of a subset, the
#' observed column is \code{NA} and a warning is given.
#'
#' @param fit A fitted model of class \code{zirt_fit}.
#'
#' @return A data frame with the columns \code{score}, \code{observed} and
#'   \code{expected} (proportions).
#'
#' @examples
#' d <- simulate_zi(n = 1500, n_items = 6, n_cat = 3, pi = 0.4, seed = 2)
#' z <- zirt_data(d$resp)
#' zi_expected(zi_fit(z))
#' @export
zi_expected <- function(fit) {
  stopifnot(inherits(fit, "zirt_fit"))
  d <- fit$data
  p <- split_par(fit$par, .sizes(fit$model, d$K))
  nodes <- if (fit$model == "hurdle") gh_grid2(fit$nodes, tanh(p$zrho)) else
    gh_nodes(fit$nodes)
  dist <- .sum_dist(.plist(p, fit$model, d$K, nodes))
  ex <- drop(crossprod(dist, nodes$w))
  if (fit$model == "zi") {
    pi0 <- stats::plogis(p$logitpi)
    ex <- (1 - pi0) * ex
    ex[1L] <- ex[1L] + pi0
  }
  complete <- stats::complete.cases(d$pat)
  if (sum(d$n[complete]) < 0.95 * sum(d$n)) {
    warning("More than 5% of respondents have missing responses; the observed ",
            "sum-score distribution is not comparable and is returned as NA.",
            call. = FALSE)
    obs <- rep(NA_real_, length(ex))
  } else {
    sc <- rowSums(d$pat[complete, , drop = FALSE])
    ob <- tapply(d$n[complete], factor(sc, levels = 0:(length(ex) - 1L)), sum)
    ob[is.na(ob)] <- 0
    obs <- as.numeric(ob) / sum(ob)
  }
  data.frame(score = 0:(length(ex) - 1L), observed = obs, expected = ex)
}
