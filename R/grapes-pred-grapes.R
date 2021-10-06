#' Binary operator for model predictions on data
#'
#' This function performs model predictions via the \code{predict} function
#' for each column of data.
#'
#' @usage data \%pred\% models
#' @param data is a matrix with rows corresponding to features, and columns
#' corresponding to cells/observations
#' @param models is a list of univariate outcome models with the features as
#' explanatory variables
#'
#' @return a matrix with rows equal to \code{length(models)} and columns
#' corresponding to cells/observations
#'
#' @keywords internal
"%pred%" <- function(data, models) {
  # data is a features x cells matrix
  # models is a list of univariate models with the features
  # as explanatory variables
  # output is a matrix with
  # rows equal to length(models) and columns are
  # cells
  do.call(rbind, lapply(models, predict, t(data)))
}
