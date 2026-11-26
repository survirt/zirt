#' Prepare item responses for zero-inflated and hurdle models
#'
#' Checks ordinal item responses, applies a scale-level filter question,
#' attaches sampling weights and design identifiers, and collapses the data
#' to distinct response patterns for fast estimation.
#'
#' @details
#' Responses must be coded 0, 1, ..., with 0 as the absence category ("no",
#' "never", "not at all"). Items may have different numbers of categories.
#' Missing responses are allowed and are ignored in the likelihood.
#'
#' A filter question is a gate asked before the scale: respondents who
#' answer "no" are not asked the items. Such respondents are outside the
#' scale by design. They are counted and reported, but they are not part of
#' the likelihood, which describes the respondents who passed the filter.
#'
#' Weights are rescaled to sum to the number of analysed respondents, so
#' that the fit is a pseudo maximum likelihood fit and information criteria
#' stay on the scale of the sample size.
#'
#' @param resp Matrix or data frame of integer item responses coded from 0.
#' @param weight Optional numeric vector of sampling weights.
#' @param filter Optional vector coded 1 (passed the filter question, items
#'   asked) and 0 (did not pass, items not asked).
#' @param strata,psu Optional stratum and primary sampling unit identifiers
#'   for design-based standard errors.
#'
#' @return An object of class \code{zirt_data}.
#'
#' @examples
#' d <- simulate_zi(n = 500, n_items = 5, n_cat = 4, pi = 0.4, seed = 1)
#' z <- zirt_data(d$resp)
#' z
#' @export
zirt_data <- function(resp, weight = NULL, filter = NULL, strata = NULL,
                      psu = NULL) {
  resp <- as.matrix(resp)
  if (!is.numeric(resp)) stop("'resp' must be numeric.", call. = FALSE)
  if (any(resp[!is.na(resp)] < 0) || any(resp[!is.na(resp)] %% 1 != 0)) {
    stop("Responses must be integers coded from 0.", call. = FALSE)
  }
  storage.mode(resp) <- "integer"
  n_all <- nrow(resp)
  if (is.null(colnames(resp))) colnames(resp) <- paste0("item", seq_len(ncol(resp)))
  if (is.null(weight)) weight <- rep(1, n_all)
  if (length(weight) != n_all) {
    stop("'weight' must have one value per row of 'resp'.", call. = FALSE)
  }
  ok <- !is.na(weight) & weight > 0
  if (!is.null(filter)) {
    if (length(filter) != n_all) {
      stop("'filter' must have one value per row of 'resp'.", call. = FALSE)
    }
    if (!all(filter[!is.na(filter)] %in% c(0, 1))) {
      stop("'filter' must be coded 0/1.", call. = FALSE)
    }
    ok <- ok & !is.na(filter)
  }
  passed <- if (is.null(filter)) rep(TRUE, n_all) else filter == 1
  answered <- rowSums(!is.na(resp)) > 0
  use <- ok & passed & answered
  if (sum(use) < 2L) stop("No analysable respondents.", call. = FALSE)
  K <- apply(resp[use, , drop = FALSE], 2L, max, na.rm = TRUE) + 1L
  if (any(K < 2L)) {
    stop("Every item needs at least two observed categories.", call. = FALSE)
  }
  r <- resp[use, , drop = FALSE]
  w <- weight[use] * sum(use) / sum(weight[use])
  key <- apply(r, 1L, paste, collapse = ",")
  ukey <- unique(key)
  idx <- match(key, ukey)
  pat <- r[match(ukey, key), , drop = FALSE]
  rownames(pat) <- NULL
  out <- list(pat = pat, n = as.numeric(rowsum(w, idx, reorder = TRUE)),
              idx = idx, w = w, K = K, items = colnames(resp),
              N = sum(use),
              zero = rowSums(pat, na.rm = TRUE) == 0L,
              strata = if (is.null(strata)) NULL else strata[use],
              psu = if (is.null(psu)) NULL else psu[use],
              filter = if (is.null(filter)) NULL else
                list(n_out = sum(ok & !passed),
                     share_out = sum(weight[ok & !passed]) / sum(weight[ok])),
              rows = which(use))
  class(out) <- "zirt_data"
  out
}

#' @export
print.zirt_data <- function(x, ...) {
  cat("Item response data for zero-inflated and hurdle models\n")
  cat(sprintf("  Respondents analysed: %d; items: %d; categories: %s\n", x$N,
              length(x$items), paste(unique(x$K), collapse = ", ")))
  cat(sprintf("  Distinct response patterns: %d\n", nrow(x$pat)))
  cat(sprintf("  All-zero respondents (weighted): %.1f%%\n",
              100 * sum(x$n[x$zero]) / sum(x$n)))
  if (!is.null(x$filter)) {
    cat(sprintf("  Stopped by the filter question: %d (weighted %.1f%% of all)\n",
                x$filter$n_out, 100 * x$filter$share_out))
  }
  if (!is.null(x$psu)) {
    cat(sprintf("  Design: %d strata, %d sampling units\n",
                if (is.null(x$strata)) 1L else length(unique(x$strata)),
                length(unique(paste(x$strata, x$psu)))))
  }
  invisible(x)
}
