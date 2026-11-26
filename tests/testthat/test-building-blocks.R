test_that("Gauss-Hermite nodes integrate normal moments exactly", {
  g <- zirt:::gh_nodes(21)
  expect_equal(sum(g$w), 1)
  expect_equal(sum(g$w * g$x), 0, tolerance = 1e-10)
  expect_equal(sum(g$w * g$x^2), 1)
  expect_equal(sum(g$w * g$x^4), 3)
  g2 <- zirt:::gh_grid2(15, 0.6)
  expect_equal(sum(g2$w * g2$t0 * g2$t1), 0.6)
  expect_equal(sum(g2$w * g2$t1^2), 1)
})

test_that("category probabilities are proper distributions", {
  th <- seq(-3, 3, length.out = 7)
  p <- zirt:::p_grm(th, 1.5, c(-0.5, 0.4, 1.2))
  expect_equal(rowSums(p), rep(1, 7))
  expect_equal(p[, 1], 1 - plogis(1.5 * (th + 0.5)))
  h <- zirt:::p_hurdle(th, th, 2, 0.3, 1.2, c(0, 1))
  expect_equal(rowSums(h), rep(1, 7))
  expect_equal(h[, 1], 1 - plogis(2 * (th - 0.3)))
  hb <- zirt:::p_hurdle(th, th, 2, 0.3, 1, numeric(0))
  expect_equal(ncol(hb), 2L)
})

test_that("the sum-score recursion matches brute-force enumeration", {
  th <- c(-1, 0.5)
  pl <- list(zirt:::p_grm(th, 1, c(0, 1)), zirt:::p_grm(th, 2, 0.5),
             zirt:::p_grm(th, 1.5, c(-0.5, 0.5)))
  d <- zirt:::.sum_dist(pl)
  expect_equal(rowSums(d), c(1, 1))
  brute <- numeric(6)
  for (a in 0:2) for (b in 0:1) for (cc in 0:2) {
    brute[a + b + cc + 1] <- brute[a + b + cc + 1] +
      pl[[1]][1, a + 1] * pl[[2]][1, b + 1] * pl[[3]][1, cc + 1]
  }
  expect_equal(d[1, ], brute)
})

test_that("zirt_data collapses patterns and handles filters and weights", {
  resp <- rbind(c(0, 0), c(0, 0), c(1, 2), c(1, 2), c(NA, 1), c(0, 1))
  z <- zirt_data(resp)
  expect_equal(nrow(z$pat), 4L)
  expect_equal(sum(z$n), 6)
  expect_equal(unname(z$K), c(2L, 3L))
  expect_equal(sum(z$zero), 1L)
  zf <- zirt_data(resp, filter = c(0, 1, 1, 1, 1, 1), weight = c(4, 1, 1, 1, 1, 1))
  expect_equal(zf$N, 5L)
  expect_equal(zf$filter$n_out, 1L)
  expect_equal(zf$filter$share_out, 4 / 9)
  expect_equal(sum(zf$n), 5)
  expect_error(zirt_data(resp - 1), "coded from 0")
  expect_error(zirt_data(resp, filter = 1:6), "0/1")
  expect_output(print(zf), "Stopped by the filter question: 1")
})

test_that("analytic gradients agree with numerical derivatives", {
  for (mod in c("grm", "zi", "hurdle")) {
    d <- simulate_zi(n = 400, n_items = 4, n_cat = 4, pi = 0.3,
                     model = if (mod == "hurdle") "hurdle" else "zi", seed = 9)
    z <- zirt_data(d$resp)
    nq <- if (mod == "hurdle") 9L else 21L
    st <- zirt:::.start(z, mod) + 0.05
    num <- zirt:::.jacobian(function(p) sum(z$n * zirt:::.ll_patterns(p, z, mod, nq)), st)
    expect_equal(zirt:::.gradient(st, z, mod, nq), drop(num), tolerance = 1e-5)
  }
})

test_that("the simulator honours its arguments", {
  d <- simulate_zi(n = 2000, n_items = 4, n_cat = 3, pi = 0.5, seed = 1)
  expect_equal(dim(d$resp), c(2000L, 4L))
  expect_true(all(rowSums(d$resp[d$structural == 1, ]) == 0))
  expect_equal(mean(d$structural), 0.5, tolerance = 0.05)
  f <- simulate_zi(n = 2000, n_items = 4, n_cat = 3, pi = 0.5,
                   filter_sens = 0.6, seed = 2)
  expect_true(all(is.na(f$resp[f$filter == 0, ])))
  expect_true(all(f$filter[f$structural == 1] == 0))
  h <- simulate_zi(n = 500, n_items = 3, n_cat = 4, model = "hurdle", seed = 3)
  expect_equal(max(h$resp), 3L)
  expect_error(simulate_zi(n_cat = 1), "at least 2")
})
