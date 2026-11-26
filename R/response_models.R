# Category probabilities -----------------------------------------------------

# Graded response model for one item: matrix nodes x K with P(Y = k | theta),
# k = 0, ..., K - 1, from a slope and K - 1 ordered thresholds.
p_grm <- function(theta, a, b) {
  cum <- stats::plogis(outer(theta, b, function(t, bb) a * (t - bb)))
  cum <- cbind(1, cum, 0)
  p <- cum[, -ncol(cum), drop = FALSE] - cum[, -1L, drop = FALSE]
  pmax(p, 1e-300)
}

# Hurdle graded response model for one item: presence follows a 2PL in the
# susceptibility trait t0, and the nonzero category follows a graded
# response model truncated at zero in the severity trait t1.
p_hurdle <- function(t0, t1, a0, b0, a1, b1) {
  pres <- stats::plogis(a0 * (t0 - b0))
  if (length(b1) == 0L) return(pmax(cbind(1 - pres, pres), 1e-300))
  sev <- p_grm(t1, a1, b1)
  pmax(cbind(1 - pres, pres * sev), 1e-300)
}

# Parameter bookkeeping ------------------------------------------------------

# Ordered thresholds from a first value and log increments.
thr_from <- function(first, logd) cumsum(c(first, exp(logd)))

# Split a parameter vector according to a named list of lengths.
split_par <- function(par, sizes) {
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  out <- lapply(seq_along(sizes), function(i)
    if (sizes[i] == 0L) numeric(0) else par[starts[i]:ends[i]])
  names(out) <- names(sizes)
  out
}

# Likelihood of every response pattern at every quadrature node, as a
# matrix patterns x nodes, from a list of per-item probability matrices.
pattern_lik <- function(pat, plist) {
  L <- matrix(1, nrow(pat), nrow(plist[[1L]]))
  for (j in seq_along(plist)) {
    y <- pat[, j]
    obs <- !is.na(y)
    if (any(obs)) {
      L[obs, ] <- L[obs, , drop = FALSE] * t(plist[[j]])[y[obs] + 1L, , drop = FALSE]
    }
  }
  L
}
