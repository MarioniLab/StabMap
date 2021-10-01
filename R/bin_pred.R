#' "%pred%"
#'
#' "%pred%"
#'
#' @param data data
#' @param models models
#'
#' @return matrix
#'
#' @examples
#'
#' @export
"%pred%" <- function(data, models) {
  # data is a features x cells matrix
  # models is a list of univariate models with the features
  # as explanatory variables
  # output is a matrix with
  # rows equal to length(models) and columns are
  # cells
  do.call(rbind, lapply(models, predict, t(data)))
}
