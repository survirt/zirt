# Model specifications -------------------------------------------------------
# Each specification gives the sizes of the blocks of the working parameter
# vector, a function returning the log-likelihood of every pattern, start
# values and the map from working to reported parameters.

.sizes <- function(model, K) {
  J <- length(K)
  if (model %in% c("grm", "zi")) {
    s <- c(loga = J, first = J, logd = sum(K - 2L))
    if (model == "zi") s <- c(s, logitpi = 1L)
  } else {
    sev <- K >= 3L
    s <- c(loga0 = J, b0 = J, loga1 = sum(sev), first1 = sum(sev),
           logd1 = sum(pmax(K - 3L, 0L)), zrho = 1L)
  }
  s
}

# Per-item probability matrices at the nodes.
.plist <- function(p, model, K, nodes) {
  J <- length(K)
  if (model %in% c("grm", "zi")) {
    off <- cumsum(c(0L, K - 2L))
    lapply(seq_len(J), function(j) {
      d <- if (K[j] > 2L) p$logd[(off[j] + 1L):off[j + 1L]] else numeric(0)
      p_grm(nodes$x, exp(p$loga[j]), thr_from(p$first[j], d))
    })
  } else {
    sev <- which(K >= 3L)
    off <- cumsum(c(0L, pmax(K[sev] - 3L, 0L)))
    lapply(seq_len(J), function(j) {
      s <- match(j, sev)
      if (is.na(s)) {
        p_hurdle(nodes$t0, nodes$t1, exp(p$loga0[j]), p$b0[j], 1, numeric(0))
      } else {
        d <- if (K[j] > 3L) p$logd1[(off[s] + 1L):off[s + 1L]] else numeric(0)
        p_hurdle(nodes$t0, nodes$t1, exp(p$loga0[j]), p$b0[j],
                 exp(p$loga1[s]), thr_from(p$first1[s], d))
      }
    })
  }
}

# Log-likelihood of every response pattern.
.ll_patterns <- function(par, d, model, nq) {
  p <- split_par(par, .sizes(model, d$K))
  if (model == "hurdle") {
    nodes <- gh_grid2(nq, tanh(p$zrho))
  } else {
    nodes <- gh_nodes(nq)
  }
  m <- drop(pattern_lik(d$pat, .plist(p, model, d$K, nodes)) %*% nodes$w)
  if (model == "zi") {
    pi0 <- stats::plogis(p$logitpi)
    m <- pi0 * d$zero + (1 - pi0) * m
  }
  log(pmax(m, 1e-300))
}

# Item to which each working parameter belongs (0 for global parameters).
.owner <- function(model, K) {
  J <- length(K)
  if (model %in% c("grm", "zi")) {
    o <- c(seq_len(J), seq_len(J), rep(seq_len(J), times = K - 2L))
    if (model == "zi") o <- c(o, 0L)
  } else {
    sev <- which(K >= 3L)
    o <- c(seq_len(J), seq_len(J), sev, sev,
           rep(sev, times = pmax(K[sev] - 3L, 0L)), 0L)
  }
  o
}

# Gradient of the weighted log-likelihood. The derivative with respect to an
# item parameter only needs the change in that item's category
# probabilities, which is obtained by a central difference at the nodes.
.gradient <- function(par, d, model, nq, eps = 1e-6) {
  K <- d$K
  J <- length(K)
  sizes <- .sizes(model, K)
  p <- split_par(par, sizes)
  nodes <- if (model == "hurdle") gh_grid2(nq, tanh(p$zrho)) else gh_nodes(nq)
  plist <- .plist(p, model, K, nodes)
  P <- nrow(d$pat)
  Q <- length(nodes$w)
  sel <- lapply(seq_len(J), function(j) {
    y <- d$pat[, j]
    S <- matrix(1, P, Q)
    obs <- !is.na(y)
    S[obs, ] <- t(plist[[j]])[y[obs] + 1L, , drop = FALSE]
    S
  })
  L <- Reduce(`*`, sel)
  m <- drop(L %*% nodes$w)
  if (model == "zi") {
    pi0 <- stats::plogis(p$logitpi)
    M <- pmax(pi0 * d$zero + (1 - pi0) * m, 1e-300)
    cf <- d$n * (1 - pi0) / M
  } else {
    cf <- d$n / pmax(m, 1e-300)
  }
  A <- sweep(L * cf, 2L, nodes$w, `*`)
  owner <- .owner(model, K)
  g <- numeric(length(par))
  for (j in seq_len(J)) {
    y <- d$pat[, j]
    obs <- !is.na(y)
    Bj <- A[obs, , drop = FALSE] / sel[[j]][obs, , drop = FALSE]
    Ck <- rowsum(Bj, factor(y[obs], levels = 0:(K[j] - 1L)), reorder = TRUE)
    full <- matrix(0, K[j], Q)
    full[as.integer(rownames(Ck)) + 1L, ] <- Ck
    for (i in which(owner == j)) {
      h <- eps * max(1, abs(par[i]))
      pp <- par; pm <- par
      pp[i] <- par[i] + h
      pm[i] <- par[i] - h
      dP <- (.plist(split_par(pp, sizes), model, K, nodes)[[j]] -
               .plist(split_par(pm, sizes), model, K, nodes)[[j]]) / (2 * h)
      g[i] <- sum(t(full) * dP)
    }
  }
  i0 <- which(owner == 0L)
  if (model == "zi") {
    g[i0] <- sum(d$n * pi0 * (1 - pi0) * (d$zero - m) / M)
  } else if (model == "hurdle") {
    h <- 1e-5
    pp <- par; pm <- par
    pp[i0] <- par[i0] + h
    pm[i0] <- par[i0] - h
    g[i0] <- (sum(d$n * .ll_patterns(pp, d, model, nq)) -
                sum(d$n * .ll_patterns(pm, d, model, nq))) / (2 * h)
  }
  g
}

# Start values from weighted category proportions.
.start <- function(d, model, pi0 = 0.2) {
  K <- d$K
  J <- length(K)
  prop <- function(j, keep) {
    y <- d$pat[, j]
    ok <- !is.na(y) & keep
    vapply(seq_len(K[j] - 1L), function(k) sum(d$n[ok & y >= k]) / sum(d$n[ok]),
           numeric(1))
  }
  clamp <- function(v) pmin(pmax(v, 0.005), 0.995)
  incr <- function(b) log(pmax(diff(b), 0.1))
  all_rows <- rep(TRUE, nrow(d$pat))
  if (model %in% c("grm", "zi")) {
    scale <- if (model == "zi") 1 - pi0 else 1
    b <- lapply(seq_len(J), function(j) -stats::qlogis(clamp(prop(j, all_rows) / scale)))
    c(rep(0, J), vapply(b, `[`, numeric(1), 1L), unlist(lapply(b, incr)),
      if (model == "zi") stats::qlogis(pi0))
  } else {
    sev <- which(K >= 3L)
    b0 <- vapply(seq_len(J), function(j) -stats::qlogis(clamp(prop(j, all_rows)[1L])),
                 numeric(1))
    b1 <- lapply(sev, function(j) {
      pos <- !is.na(d$pat[, j]) & d$pat[, j] >= 1L
      cp <- vapply(2:(K[j] - 1L), function(k)
        sum(d$n[pos & d$pat[, j] >= k]) / sum(d$n[pos]), numeric(1))
      -stats::qlogis(clamp(cp))
    })
    c(rep(0, J), b0, rep(0, length(sev)), vapply(b1, `[`, numeric(1), 1L),
      unlist(lapply(b1, incr)), atanh(0.5))
  }
}

# Reported parameters on their natural scale, as a named vector.
.natural <- function(par, d, model) {
  K <- d$K
  J <- length(K)
  p <- split_par(par, .sizes(model, K))
  out <- numeric(0)
  if (model %in% c("grm", "zi")) {
    off <- cumsum(c(0L, K - 2L))
    for (j in seq_len(J)) {
      dd <- if (K[j] > 2L) p$logd[(off[j] + 1L):off[j + 1L]] else numeric(0)
      b <- thr_from(p$first[j], dd)
      v <- c(exp(p$loga[j]), b)
      names(v) <- paste0(d$items[j], "|", c("a", paste0("b", seq_along(b))))
      out <- c(out, v)
    }
    if (model == "zi") out <- c(out, "structural zeros|pi" = stats::plogis(p$logitpi))
  } else {
    sev <- which(K >= 3L)
    off <- cumsum(c(0L, pmax(K[sev] - 3L, 0L)))
    for (j in seq_len(J)) {
      v <- c(exp(p$loga0[j]), p$b0[j])
      names(v) <- paste0(d$items[j], "|", c("a0", "b0"))
      s <- match(j, sev)
      if (!is.na(s)) {
        dd <- if (K[j] > 3L) p$logd1[(off[s] + 1L):off[s + 1L]] else numeric(0)
        b <- thr_from(p$first1[s], dd)
        v1 <- c(exp(p$loga1[s]), b)
        names(v1) <- paste0(d$items[j], "|", c("a1", paste0("b1_", seq_along(b))))
        v <- c(v, v1)
      }
      out <- c(out, v)
    }
    out <- c(out, "susceptibility-severity|rho" = tanh(p$zrho))
  }
  out
}

# Central-difference Jacobian of a vector function.
.jacobian <- function(f, x, eps = 1e-5) {
  f0 <- f(x)
  J <- matrix(0, length(f0), length(x))
  for (i in seq_along(x)) {
    h <- eps * max(1, abs(x[i]))
    xp <- x; xm <- x
    xp[i] <- x[i] + h
    xm[i] <- x[i] - h
    J[, i] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}

# Covariance of the working parameters: model-based (inverse Hessian) or
# design-based (Taylor linearisation with strata and sampling units).
.vcov <- function(par, d, model, nq, hess, type) {
  Hinv <- tryCatch(solve(hess), error = function(e) NULL)
  if (is.null(Hinv)) return(matrix(NA_real_, length(par), length(par)))
  if (type == "model") return(Hinv)
  S <- .jacobian(function(p) .ll_patterns(p, d, model, nq), par)
  sc <- S[d$idx, , drop = FALSE] * d$w
  strata <- if (is.null(d$strata)) rep(1L, d$N) else d$strata
  psu <- if (is.null(d$psu)) seq_len(d$N) else d$psu
  cl <- paste(strata, psu, sep = "\r")
  z <- rowsum(sc, cl, reorder = FALSE)
  zs <- strata[match(rownames(z), cl)]
  G <- matrix(0, length(par), length(par))
  for (h in unique(zs)) {
    zh <- z[zs == h, , drop = FALSE]
    nh <- nrow(zh)
    if (nh < 2L) next
    zc <- sweep(zh, 2L, colMeans(zh))
    G <- G + nh / (nh - 1) * crossprod(zc)
  }
  Hinv %*% G %*% Hinv
}

# Shared estimation routine.
.fit <- function(d, model, nq, starts, vcov, control) {
  stopifnot(inherits(d, "zirt_data"))
  if (model == "hurdle" && !any(d$K >= 3L)) {
    stop("The hurdle model needs at least one item with three or more categories.",
         call. = FALSE)
  }
  nll <- function(par) -sum(d$n * .ll_patterns(par, d, model, nq))
  ngr <- function(par) -.gradient(par, d, model, nq)
  inits <- if (model == "zi") {
    p0 <- sum(d$n[d$zero]) / sum(d$n)
    lapply(seq_len(starts), function(s) .start(d, "zi", pi0 = p0 * c(0.5, 0.2, 0.8, 0.05, 0.95)[(s - 1L) %% 5L + 1L]))
  } else list(.start(d, model))
  ctrl <- utils::modifyList(list(maxit = 2000L, reltol = 1e-10), control)
  best <- NULL
  for (st in inits) {
    o <- tryCatch(stats::optim(st, nll, ngr, method = "BFGS", control = ctrl),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) stop("Estimation failed for all start values.", call. = FALSE)
  hess <- stats::optimHess(best$par, nll, ngr)
  V <- .vcov(best$par, d, model, nq, hess, vcov)
  nat <- .natural(best$par, d, model)
  Jn <- .jacobian(function(p) .natural(p, d, model), best$par)
  se <- sqrt(pmax(diag(Jn %*% V %*% t(Jn)), 0))
  lab <- strsplit(names(nat), "|", fixed = TRUE)
  coef <- data.frame(item = vapply(lab, `[`, "", 1L),
                     parameter = vapply(lab, `[`, "", 2L),
                     estimate = unname(nat), se = se, stringsAsFactors = FALSE)
  npar <- length(best$par)
  out <- list(model = model, par = best$par, coefficients = coef,
              logLik = -best$value, npar = npar, n = d$N,
              AIC = 2 * best$value + 2 * npar,
              BIC = 2 * best$value + log(d$N) * npar,
              convergence = best$convergence, vcov = V, vcov_type = vcov,
              nodes = nq, data = d)
  class(out) <- "zirt_fit"
  out
}
