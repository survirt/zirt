#' Fit the graded response model
#'
#' Marginal maximum likelihood estimation of Samejima's graded response
#' model (the two-parameter logistic model for binary items) with a standard
#' normal latent trait. It is the reference model against which the
#' zero-inflated and hurdle models are compared.
#'
#' @param data An object created by \code{\link{zirt_data}}.
#' @param nodes Number of Gauss-Hermite quadrature nodes.
#' @param vcov \code{"model"} for standard errors from the inverse Hessian,
#'   or \code{"design"} for Taylor-linearisation standard errors that use the
#'   weights, strata and sampling units stored in \code{data}.
#' @param control List of control parameters passed to
#'   \code{\link[stats]{optim}}.
#'
#' @return An object of class \code{zirt_fit}.
#'
#' @references Samejima, F. (1969). Estimation of latent ability using a
#'   response pattern of graded scores. \emph{Psychometrika Monograph
#'   Supplement, 34}(4, Pt. 2), 1-97.
#'
#' @examples
#' d <- simulate_zi(n = 800, n_items = 5, n_cat = 3, pi = 0, seed = 1)
#' fit <- grm_fit(zirt_data(d$resp))
#' fit
#' @export
grm_fit <- function(data, nodes = 41L, vcov = c("model", "design"),
                    control = list()) {
  .fit(data, "grm", as.integer(nodes), 1L, match.arg(vcov), control)
}

#' Fit the zero-inflated mixture graded response model
#'
#' Estimates a two-class mixture in which a proportion \eqn{\pi} of
#' respondents are structural zeros, who answer in the absence category on
#' every item with probability one, and the remaining respondents follow a
#' graded response model with a standard normal latent trait.
#'
#' @details
#' The marginal probability of response pattern \eqn{y} is
#' \deqn{P(y) = \pi\, 1\{y = 0\} + (1 - \pi) \int \prod_j P_j(y_j \mid
#' \theta)\, \phi(\theta)\, d\theta .}
#' The latent trait is scaled in the class that has the trait, so item
#' parameters describe respondents with some level of the construct (Wall,
#' Park, & Moustaki, 2015). With binary items the model is the zero-inflated
#' two-parameter logistic model. Because mixture likelihoods can have local
#' maxima, estimation is repeated from several start values of \eqn{\pi} and
#' the solution with the highest likelihood is returned.
#'
#' @inheritParams grm_fit
#' @param starts Number of start values for the proportion of structural
#'   zeros.
#'
#' @return An object of class \code{zirt_fit}.
#'
#' @references Wall, M. M., Park, J. Y., & Moustaki, I. (2015). IRT modeling
#'   in the presence of zero-inflation with application to psychiatric
#'   disorder severity. \emph{Applied Psychological Measurement, 39}(8),
#'   583-597.
#'
#' @examples
#' d <- simulate_zi(n = 1500, n_items = 6, n_cat = 3, pi = 0.4, seed = 2)
#' fit <- zi_fit(zirt_data(d$resp))
#' fit
#' @export
zi_fit <- function(data, nodes = 41L, starts = 3L,
                   vcov = c("model", "design"), control = list()) {
  .fit(data, "zi", as.integer(nodes), as.integer(starts), match.arg(vcov),
       control)
}

#' Fit the multidimensional hurdle graded response model
#'
#' Estimates a two-part model in which each item has a presence component
#' and a severity component that depend on two correlated latent traits.
#'
#' @details
#' Whether the absence category is left follows a two-parameter logistic
#' model in the susceptibility trait \eqn{\theta_0}. Given a nonzero
#' response, the category follows a graded response model truncated at zero
#' in the severity trait \eqn{\theta_1}:
#' \deqn{P(Y_j = 0) = 1 - \Psi\{a_{0j}(\theta_0 - b_{0j})\}, \qquad
#' P(Y_j = k) = \Psi\{a_{0j}(\theta_0 - b_{0j})\}\, T_j(k \mid \theta_1),
#' \quad k \ge 1 .}
#' The two traits are bivariate standard normal with correlation
#' \eqn{\rho} (Magnus & Liu, 2022; Magnus, 2024). The likelihood is
#' integrated on a two-dimensional Gauss-Hermite grid. Binary items
#' contribute to the presence component only. At least one item needs three
#' or more categories.
#'
#' @inheritParams grm_fit
#' @param nodes Number of quadrature nodes per dimension.
#'
#' @return An object of class \code{zirt_fit}.
#'
#' @references
#' Magnus, B. E. (2024). Item response modeling of clinical instruments with
#' filter questions: Disentangling symptom presence and severity.
#' \emph{Applied Psychological Measurement, 48}(6), 235-256.
#'
#' Magnus, B. E., & Liu, Y. (2022). Symptom presence and symptom severity as
#' unique indicators of psychopathology: An application of multidimensional
#' zero-inflated and hurdle graded response models. \emph{Educational and
#' Psychological Measurement, 82}(5), 938-966.
#'
#' @examples
#' \donttest{
#' d <- simulate_zi(n = 1500, n_items = 5, n_cat = 4, model = "hurdle",
#'                  rho = 0.5, seed = 3)
#' fit <- hurdle_fit(zirt_data(d$resp), nodes = 15)
#' fit
#' }
#' @export
hurdle_fit <- function(data, nodes = 21L, vcov = c("model", "design"),
                       control = list()) {
  .fit(data, "hurdle", as.integer(nodes), 1L, match.arg(vcov), control)
}

#' Compare zero-inflated, hurdle and standard models
#'
#' Tabulates the log-likelihood, number of parameters, AIC and BIC of fitted
#' models and tests the graded response model against the zero-inflated
#' model.
#'
#' @details
#' The graded response model is the zero-inflated model with \eqn{\pi = 0},
#' a value on the boundary of the parameter space. The likelihood-ratio
#' statistic is therefore referred to an equal mixture of a point mass at
#' zero and a chi-square distribution with one degree of freedom (Self &
#' Liang, 1987), and the p value is half the usual one. The hurdle model is
#' not nested in the others and is compared with information criteria. With
#' sampling weights the fits are pseudo likelihoods and the test is
#' approximate.
#'
#' @param ... Objects of class \code{zirt_fit} fitted to the same data.
#'
#' @return A data frame of class \code{zirt_compare}.
#'
#' @references Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of
#'   maximum likelihood estimators and likelihood ratio tests under
#'   nonstandard conditions. \emph{Journal of the American Statistical
#'   Association, 82}(398), 605-610.
#'
#' @examples
#' d <- simulate_zi(n = 1500, n_items = 6, n_cat = 3, pi = 0.4, seed = 2)
#' z <- zirt_data(d$resp)
#' zi_compare(grm_fit(z), zi_fit(z))
#' @export
zi_compare <- function(...) {
  fits <- list(...)
  if (!all(vapply(fits, inherits, logical(1), "zirt_fit"))) {
    stop("All arguments must be fitted zirt models.", call. = FALSE)
  }
  lab <- c(grm = "Graded response", zi = "Zero-inflated mixture",
           hurdle = "Multidimensional hurdle")
  tab <- data.frame(model = unname(lab[vapply(fits, `[[`, "", "model")]),
                    logLik = vapply(fits, `[[`, 0, "logLik"),
                    npar = vapply(fits, `[[`, 0L, "npar"),
                    AIC = vapply(fits, `[[`, 0, "AIC"),
                    BIC = vapply(fits, `[[`, 0, "BIC"),
                    stringsAsFactors = FALSE)
  tab$LR <- NA_real_
  tab$p <- NA_real_
  m <- vapply(fits, `[[`, "", "model")
  if (all(c("grm", "zi") %in% m)) {
    lr <- max(0, 2 * (tab$logLik[m == "zi"][1L] - tab$logLik[m == "grm"][1L]))
    tab$LR[m == "zi"] <- lr
    tab$p[m == "zi"] <- 0.5 * stats::pchisq(lr, 1L, lower.tail = FALSE)
  }
  class(tab) <- c("zirt_compare", "data.frame")
  tab
}
