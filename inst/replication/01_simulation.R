# Simulation study for zirt.
# Part A: data from the zero-inflated mixture model; fit the graded response
#         and zero-inflated models.
# Part B: data from the hurdle model; fit all three models.
#
# Usage: Rscript 01_simulation.R <output dir> [reps A] [reps B] [cores]

library(zirt)
library(parallel)

args <- commandArgs(TRUE)
out_dir <- args[1]
RA <- if (length(args) > 1) as.integer(args[2]) else 200L
RB <- if (length(args) > 2) as.integer(args[3]) else 100L
cores <- if (length(args) > 3) as.integer(args[4]) else max(1L, detectCores() - 1L)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

condA <- expand.grid(pi = c(0, 0.3, 0.6), form = c("4 binary", "6 four-category"),
                     n = c(1000L, 5000L), stringsAsFactors = FALSE)
condA$cond <- seq_len(nrow(condA))

repA <- function(job) {
  cd <- job$cd
  set.seed(100000L * cd$cond + job$r)
  binary <- cd$form == "4 binary"
  d <- if (binary) simulate_zi(n = cd$n, n_items = 4, n_cat = 2, pi = cd$pi,
                               a = c(1.5, 2, 2.5, 1.8), b = list(0.5, 0.8, 1.0, 0.3)) else
    simulate_zi(n = cd$n, n_items = 6, n_cat = 4, pi = cd$pi)
  z <- zirt_data(d$resp)
  g <- grm_fit(z, nodes = 31)
  f <- zi_fit(z, nodes = 31, starts = 2)
  cf <- coef(f); cg <- coef(g)
  pr <- cf[cf$parameter == "pi", ]
  c(cond = cd$cond, rep = job$r, pi_hat = pr$estimate, pi_se = pr$se,
    cover = as.numeric(abs(pr$estimate - cd$pi) < qnorm(.975) * pr$se),
    a_bias_zi = mean(cf$estimate[cf$parameter == "a"] / d$a) - 1,
    a_bias_grm = mean(cg$estimate[cg$parameter == "a"] / d$a) - 1,
    b_bias_zi = mean(cf$estimate[cf$parameter == "b1"] - vapply(d$b, `[`, 0, 1)),
    b_bias_grm = mean(cg$estimate[cg$parameter == "b1"] - vapply(d$b, `[`, 0, 1)),
    lr_reject = as.numeric(zi_compare(g, f)$p[2] < .05),
    bic_zi = as.numeric(f$BIC < g$BIC),
    conv = as.numeric(f$convergence == 0 && g$convergence == 0))
}

condB <- data.frame(rho = c(0.3, 0.7), cond = 1:2)
repB <- function(job) {
  cd <- job$cd
  set.seed(900000L + 100000L * cd$cond + job$r)
  d <- simulate_zi(n = 2000, n_items = 5, n_cat = 4, model = "hurdle", rho = cd$rho)
  z <- zirt_data(d$resp)
  g <- grm_fit(z, nodes = 31)
  f <- zi_fit(z, nodes = 31, starts = 2)
  h <- hurdle_fit(z, nodes = 11)
  ch <- coef(h)
  rr <- ch[ch$parameter == "rho", ]
  bic <- c(g$BIC, f$BIC, h$BIC)
  c(cond = cd$cond, rep = job$r, rho_hat = rr$estimate, rho_se = rr$se,
    cover = as.numeric(abs(rr$estimate - cd$rho) < qnorm(.975) * rr$se),
    b0_bias = mean(ch$estimate[ch$parameter == "b0"] - vapply(d$b, `[`, 0, 1)),
    a0_bias = mean(ch$estimate[ch$parameter == "a0"] / d$a) - 1,
    best_hurdle = as.numeric(which.min(bic) == 3L),
    best_zi = as.numeric(which.min(bic) == 2L),
    best_grm = as.numeric(which.min(bic) == 1L))
}

mk <- function(conds, R) unlist(lapply(seq_len(nrow(conds)), function(i)
  lapply(seq_len(R), function(r) list(cd = conds[i, ], r = r))), recursive = FALSE)
t0 <- Sys.time()
cl <- makeCluster(cores)
clusterEvalQ(cl, library(zirt))
resA <- as.data.frame(do.call(rbind, parLapplyLB(cl, mk(condA, RA), repA)))
tA <- Sys.time()
resB <- as.data.frame(do.call(rbind, parLapplyLB(cl, mk(condB, RB), repB)))
stopCluster(cl)

write.csv(resA, file.path(out_dir, "simulation_A_replications.csv"), row.names = FALSE)
write.csv(resB, file.path(out_dir, "simulation_B_replications.csv"), row.names = FALSE)
sumA <- merge(condA, aggregate(. ~ cond, data = resA[, names(resA) != "rep"], FUN = mean), by = "cond")
sumA$pi_rmse <- sqrt(tapply((resA$pi_hat - condA$pi[resA$cond])^2, resA$cond, mean))[as.character(sumA$cond)]
sumA$replications <- RA
sumB <- merge(condB, aggregate(. ~ cond, data = resB[, names(resB) != "rep"], FUN = mean), by = "cond")
sumB$replications <- RB
write.csv(sumA, file.path(out_dir, "simulation_A_summary.csv"), row.names = FALSE)
write.csv(sumB, file.path(out_dir, "simulation_B_summary.csv"), row.names = FALSE)
cat("Part A:", RA, "replications x", nrow(condA), "conditions,",
    round(as.numeric(difftime(tA, t0, units = "mins")), 1), "min\n")
cat("Part B:", RB, "replications x", nrow(condB), "conditions,",
    round(as.numeric(difftime(Sys.time(), tA, units = "mins")), 1), "min\n")
print(round(sumA[, c("pi", "n", "pi_hat", "pi_rmse", "cover", "a_bias_zi", "a_bias_grm", "lr_reject", "bic_zi")], 3))
print(round(sumB[, -1], 3))
sink(file.path(out_dir, "sessionInfo_simulation.txt")); print(sessionInfo()); sink()
