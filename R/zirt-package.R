#' zirt: Zero-inflated and hurdle item response models
#'
#' Item response models for scales on which most respondents report no
#' symptoms or problems, with support for filter questions and survey
#' designs.
#'
#' Prepare the data with \code{\link{zirt_data}}; fit
#' \code{\link{grm_fit}}, \code{\link{zi_fit}} and
#' \code{\link{hurdle_fit}}; choose a model with \code{\link{zi_compare}}
#' and \code{\link{zi_expected}}; and obtain scores with
#' \code{\link{zi_score}}. \code{\link{simulate_zi}} generates data from
#' the zero-inflated and hurdle models.
#'
#' @keywords internal
#' @importFrom utils modifyList
"_PACKAGE"
