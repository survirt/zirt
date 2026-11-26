test_that("grm_fit reproduces mirt", {
  skip_if_not_installed("mirt")
  d <- simulate_zi(n = 1500, n_items = 5, n_cat = 3, pi = 0, seed = 11)
  fit <- grm_fit(zirt_data(d$resp))
  ref <- mirt::mirt(as.data.frame(d$resp), 1, itemtype = "graded",
                    verbose = FALSE, quadpts = 41, TOL = 1e-6)
  expect_equal(fit$logLik, mirt::extract.mirt(ref, "logLik"), tolerance = 1e-5)
  cm <- mirt::coef(ref, IRTpars = TRUE, simplify = TRUE)$items
  est <- matrix(coef(fit)$estimate, 5, 3, byrow = TRUE)
  expect_lt(max(abs(est - cm)), 0.01)
})

test_that("weights act as frequencies", {
  d <- simulate_zi(n = 300, n_items = 4, n_cat = 2, pi = 0, seed = 12)
  w <- rep(c(1, 2), 150)
  a <- grm_fit(zirt_data(d$resp, weight = w), nodes = 21)
  dup <- d$resp[rep(seq_len(300), times = w), ]
  b <- grm_fit(zirt_data(dup), nodes = 21)
  expect_equal(coef(a)$estimate, coef(b)$estimate, tolerance = 1e-4)
  expect_equal(a$logLik * 450 / 300, b$logLik, tolerance = 1e-4)
})

test_that("zi_fit recovers the proportion of structural zeros", {
  d <- simulate_zi(n = 3000, n_items = 6, n_cat = 3, pi = 0.4, seed = 13)
  z <- zirt_data(d$resp)
  g <- grm_fit(z, nodes = 31)
  f <- zi_fit(z, nodes = 31, starts = 2)
  cf <- coef(f)
  est <- cf$estimate[cf$parameter == "pi"]
  expect_lt(abs(est - mean(d$structural)), 0.06)
  expect_gt(f$logLik, g$logLik)
  cmp <- zi_compare(g, f)
  expect_s3_class(cmp, "zirt_compare")
  expect_lt(cmp$p[2], 0.001)
  expect_equal(cmp$model[which.min(cmp$BIC)], "Zero-inflated mixture")
  # the standard model overstates the slopes when zeros are inflated
  expect_gt(mean(coef(g)$estimate[coef(g)$parameter == "a"]),
            mean(cf$estimate[cf$parameter == "a"]))
  sc <- zi_score(f)
  expect_equal(nrow(sc), 3000L)
  expect_true(all(sc$p_structural >= 0 & sc$p_structural <= 1))
  expect_true(all(sc$p_structural[rowSums(d$resp) > 0] == 0))
  expect_equal(mean(sc$p_structural), est, tolerance = 1e-3)
  ex <- zi_expected(f)
  expect_equal(sum(ex$expected), 1)
  expect_equal(ex$expected[1], ex$observed[1], tolerance = 0.01)
})

test_that("without inflation the zero-inflated model returns to the standard model", {
  d <- simulate_zi(n = 1500, n_items = 5, n_cat = 3, pi = 0, seed = 14)
  z <- zirt_data(d$resp)
  g <- grm_fit(z, nodes = 31)
  f <- zi_fit(z, nodes = 31, starts = 2)
  expect_lt(f$logLik - g$logLik, 2)
  expect_lt(coef(f)$estimate[coef(f)$parameter == "pi"], 0.1)
  expect_gt(zi_compare(g, f)$p[2], 0.01)
})

test_that("hurdle_fit separates susceptibility and severity", {
  skip_on_cran()
  d <- simulate_zi(n = 2000, n_items = 4, n_cat = 4, model = "hurdle",
                   rho = 0.5, seed = 15)
  z <- zirt_data(d$resp)
  h <- hurdle_fit(z, nodes = 11)
  g <- grm_fit(z, nodes = 31)
  cf <- coef(h)
  expect_lt(abs(cf$estimate[cf$parameter == "rho"] - 0.5), 0.2)
  expect_lt(max(abs(cf$estimate[cf$parameter == "b0"] - vapply(d$b, `[`, 0, 1))), 0.3)
  expect_lt(h$BIC, g$BIC)
  sc <- zi_score(h, patterns = TRUE)
  expect_true(all(c("susceptibility", "severity", "n") %in% names(sc)))
  expect_equal(sum(sc$n), 2000)
  expect_error(hurdle_fit(zirt_data(d$resp > 0)), "three or more")
})

test_that("design-based standard errors use strata and sampling units", {
  d <- simulate_zi(n = 1200, n_items = 4, n_cat = 3, pi = 0.3, seed = 16)
  psu <- rep(1:60, each = 20)
  strata <- rep(1:30, each = 40)
  z <- zirt_data(d$resp, weight = runif(1200, 0.5, 2), strata = strata, psu = psu)
  a <- zi_fit(z, nodes = 21, starts = 1)
  b <- zi_fit(z, nodes = 21, starts = 1, vcov = "design")
  expect_equal(coef(a)$estimate, coef(b)$estimate)
  expect_true(all(is.finite(coef(b)$se)))
  expect_true(all(eigen(b$vcov, symmetric = TRUE)$values > -1e-8))
  expect_equal(b$vcov_type, "design")
})

test_that("methods print and plot without error", {
  d <- simulate_zi(n = 800, n_items = 4, n_cat = 3, pi = 0.3,
                   filter_sens = 0.7, seed = 17)
  f <- zi_fit(zirt_data(d$resp, filter = d$filter), nodes = 21, starts = 1)
  expect_output(print(f), "Zero-inflated mixture")
  expect_output(print(summary(f)), "filter question")
  expect_s3_class(logLik(f), "logLik")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(f))
  expect_error(zi_compare(f, 1), "fitted zirt models")
})

test_that("zi_expected flags data with planned missingness", {
  d <- simulate_zi(n = 800, n_items = 4, n_cat = 3, pi = 0.3, seed = 18)
  resp <- d$resp
  resp[rowSums(resp[, 1:2]) == 0, 3:4] <- NA   # follow-up items not asked
  f <- zi_fit(zirt_data(resp), nodes = 21, starts = 1)
  expect_warning(e <- zi_expected(f), "missing responses")
  expect_true(all(is.na(e$observed)))
  expect_equal(sum(e$expected), 1)
  full <- zi_expected(zi_fit(zirt_data(d$resp), nodes = 21, starts = 1))
  expect_false(anyNA(full$observed))
})
