# zirt

Item response models for scales where most people score zero.

Symptom checklists, adherence items and hardship scales in population health
surveys share two features. A large part of the sample reports nothing at
all, and many scales sit behind a filter question, so that only some
respondents are asked the items. A standard item response model treats every
zero as a low trait level and lets the pile of zeros distort slopes,
thresholds and scores. `zirt` fits three models side by side and lets the
data choose:

- the graded response model (`grm_fit`);
- the zero-inflated mixture graded response model (`zi_fit`), which adds a
  class of structural zeros (Wall, Park, & Moustaki, 2015);
- the multidimensional hurdle graded response model (`hurdle_fit`), which
  separates the susceptibility to have a symptom from its severity (Magnus &
  Liu, 2022; Magnus, 2024).

Binary and ordinal items, missing responses, scale-level filter questions,
sampling weights and stratified cluster designs are supported. Estimation is
marginal maximum likelihood with Gauss-Hermite quadrature and analytic
gradients, written in base R.

## Install

    remotes::install_github("survirt/zirt")

R 4.1.0 or later on any platform. No compiled code; imports `stats`,
`graphics` and `utils` only.

## A short session

    library(zirt)

    sim <- simulate_zi(n = 3000, n_items = 6, n_cat = 4, pi = 0.4, seed = 1)
    dat <- zirt_data(sim$resp)

    m0 <- grm_fit(dat)
    m1 <- zi_fit(dat)
    m2 <- hurdle_fit(dat)

    zi_compare(m0, m1, m2)     # log-likelihood, AIC, BIC, boundary LR test
    summary(m1)                # item parameters, share of structural zeros
    plot(m1)                   # observed and implied sum-score distributions
    head(zi_score(m1))         # trait score and P(structural zero)

With a filter question and a survey design:

    dat <- zirt_data(resp, weight = w, filter = took_medication,
                     strata = stratum, psu = cluster)
    fit <- zi_fit(dat, vcov = "design")

## What the functions do

    zirt_data()     check, weight and collapse responses; apply a filter
    grm_fit()       graded response model (2PL for binary items)
    zi_fit()        zero-inflated mixture model
    hurdle_fit()    multidimensional hurdle model
    zi_compare()    fit table and likelihood-ratio test on the boundary
    zi_score()      EAP scores, P(structural zero), susceptibility, severity
    zi_expected()   observed and model-implied sum-score distributions
    simulate_zi()   data generator, with an optional filter question

Scripts for the simulation study and for the analysis of the 2023 National
Health Interview Survey are in `inst/replication`. The survey files are public
and are downloaded from the National Center for Health Statistics.

MIT license.
