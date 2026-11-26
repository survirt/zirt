# Parametric bootstrap at the fitted values of the cost-related medication
# nonadherence scale.
#
# Usage: Rscript 03_bootstrap_medication.R <output dir> [replications]

library(zirt)

args <- commandArgs(TRUE)
out_dir <- args[1]
R <- if (length(args) > 1) as.integer(args[2]) else 200L
cf <- read.csv(file.path(out_dir, "nhis_coefficients.csv"), stringsAsFactors = FALSE)
cf <- cf[cf$scale == "medication" & cf$model == "zi", ]
a <- cf$estimate[cf$parameter == "a"]
b <- as.list(cf$estimate[cf$parameter == "b1"])
pi0 <- cf$estimate[cf$parameter == "pi"]
n <- read.csv(file.path(out_dir, "nhis_zeros.csv"))
n <- n$n[n$scale == "medication"]

res <- t(vapply(seq_len(R), function(r) {
  set.seed(770000L + r)
  d <- simulate_zi(n = n, n_items = 4, n_cat = 2, pi = pi0, a = a, b = b)
  z <- zirt_data(d$resp)
  g <- grm_fit(z)
  f <- zi_fit(z, starts = 5)
  pr <- coef(f)[coef(f)$parameter == "pi", ]
  c(rep = r, pi_hat = pr$estimate, pi_se = pr$se,
    cover = as.numeric(abs(pr$estimate - pi0) < qnorm(.975) * pr$se),
    lr_reject = as.numeric(zi_compare(g, f)$p[2] < .05),
    a_ratio_zi = mean(coef(f)$estimate[coef(f)$parameter == "a"] / a),
    a_ratio_grm = mean(coef(g)$estimate[coef(g)$parameter == "a"] / a),
    all_zero = mean(rowSums(d$resp) == 0))
}, numeric(8)))
res <- as.data.frame(res)
write.csv(res, file.path(out_dir, "bootstrap_medication_replications.csv"), row.names = FALSE)
s <- data.frame(replications = R, n = n, pi_true = pi0, pi_mean = mean(res$pi_hat),
                pi_sd = sd(res$pi_hat), pi_se_mean = mean(res$pi_se, na.rm = TRUE),
                pi_q05 = unname(quantile(res$pi_hat, .05)), pi_q95 = unname(quantile(res$pi_hat, .95)),
                cover = mean(res$cover, na.rm = TRUE), lr_reject = mean(res$lr_reject),
                a_ratio_zi = mean(res$a_ratio_zi), a_ratio_grm = mean(res$a_ratio_grm),
                all_zero = mean(res$all_zero))
write.csv(s, file.path(out_dir, "bootstrap_medication_summary.csv"), row.names = FALSE)
print(round(s, 4))
