# National Health Interview Survey 2023, sample adult file: three scales
# with excess zeros and filter questions.
#
# Usage: Rscript 02_nhis2023.R <adult23.csv> <output dir>
#
# The public-use file adult23.csv is available from the National Center for
# Health Statistics:
# https://www.cdc.gov/nchs/nhis/documentation/2023-nhis.html

library(zirt)

args <- commandArgs(TRUE)
d <- read.csv(args[1])
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

yn <- function(v) ifelse(v == 1, 1L, ifelse(v == 2, 0L, NA_integer_))
w <- d$WTFA_A

scales <- list()

# 1. PHQ-4: depressive and anxiety symptoms, past two weeks
#    1 not at all ... 4 nearly every day -> 0 ... 3
phq <- sapply(c("PHQ41_A", "PHQ42_A", "PHQ43_A", "PHQ44_A"), function(v)
  ifelse(d[[v]] %in% 1:4, d[[v]] - 1L, NA_integer_))
colnames(phq) <- c("interest", "depressed", "anxious", "worry")
scales$phq4 <- zirt_data(phq, weight = w, strata = d$PSTRAT, psu = d$PPSU)

# 2. Cost-related medication nonadherence, past 12 months. Three items are
#    asked only of adults who took prescription medication (filter question).
rx <- cbind(skipped_doses = yn(d$RXSK12M_A), took_less = yn(d$RXLS12M_A),
            delayed_filling = yn(d$RXDL12M_A), did_not_get = yn(d$RXDG12M_A))
scales$medication <- zirt_data(rx, weight = w, filter = yn(d$RX12M_A),
                               strata = d$PSTRAT, psu = d$PPSU)

# 3. Household food insecurity, past 30 days. Three screening statements
#    (never, sometimes, often true) and five follow-up items asked when a
#    screening statement was affirmed.
tri <- function(v) ifelse(v %in% 1:3, 3L - v, NA_integer_)
fds <- cbind(worried = tri(d$FDSRUNOUT_A), did_not_last = tri(d$FDSLAST_A),
             unbalanced = tri(d$FDSBALANCE_A), cut_meals = yn(d$FDSSKIP_A),
             ate_less = yn(d$FDSLESS_A), hungry = yn(d$FDSHUNGRY_A),
             lost_weight = yn(d$FDSWEIGHT_A), no_food_day = yn(d$FDSNOTEAT_A))
scales$food <- zirt_data(fds, weight = w, strata = d$PSTRAT, psu = d$PPSU)

timing <- list(); comp <- list(); coefs <- list(); zeros <- list()
expd <- list(); scor <- list()
for (s in names(scales)) {
  z <- scales[[s]]
  t0 <- Sys.time()
  fits <- list(grm = grm_fit(z, vcov = "design"),
               zi = zi_fit(z, starts = 5, vcov = "design"))
  if (any(z$K >= 3L)) fits$hurdle <- hurdle_fit(z, nodes = 21, vcov = "design")
  timing[[s]] <- data.frame(scale = s, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cm <- do.call(zi_compare, unname(fits))
  comp[[s]] <- cbind(scale = s, n = z$N, patterns = nrow(z$pat), as.data.frame(cm),
                     converged = vapply(fits, function(f) f$convergence == 0L, logical(1)))
  coefs[[s]] <- do.call(rbind, lapply(names(fits), function(m)
    cbind(scale = s, model = m, coef(fits[[m]]))))
  pi_row <- coef(fits$zi)[coef(fits$zi)$parameter == "pi", ]
  zeros[[s]] <- data.frame(scale = s, n = z$N,
    filtered_out = if (is.null(z$filter)) NA else z$filter$share_out,
    all_zero = sum(z$n[z$zero]) / sum(z$n), pi = pi_row$estimate, pi_se = pi_row$se)
  expd[[s]] <- do.call(rbind, lapply(names(fits), function(m)
    cbind(scale = s, model = m, zi_expected(fits[[m]]))))
  zp <- which(z$zero)[1]
  sg <- zi_score(fits$grm, patterns = TRUE)
  sz <- zi_score(fits$zi, patterns = TRUE)
  scor[[s]] <- data.frame(scale = s, share_zero_pattern = z$n[zp] / sum(z$n),
    theta_zero_grm = sg$theta[zp], sd_zero_grm = sg$sd[zp],
    theta_zero_zi = sz$theta[zp], p_structural_zero_pattern = sz$p_structural[zp],
    cor_scores = stats::cov.wt(cbind(sg$theta, sz$theta), wt = z$n / sum(z$n), cor = TRUE)$cor[1, 2])
  message(s, ": done")
}
write.csv(do.call(rbind, comp), file.path(out_dir, "nhis_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, coefs), file.path(out_dir, "nhis_coefficients.csv"), row.names = FALSE)
write.csv(do.call(rbind, zeros), file.path(out_dir, "nhis_zeros.csv"), row.names = FALSE)
write.csv(do.call(rbind, expd), file.path(out_dir, "nhis_expected.csv"), row.names = FALSE)
write.csv(do.call(rbind, scor), file.path(out_dir, "nhis_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, timing), file.path(out_dir, "nhis_timing.csv"), row.names = FALSE)
sink(file.path(out_dir, "sessionInfo_nhis.txt")); print(sessionInfo()); sink()
