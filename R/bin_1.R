#' "%*1%"
#'
#' "%*1%"
#'
#' @param a a
#' @param b b
#'
#' @return matrix
#'
#' @examples
#'
#' @export
"%*1%" <- function(a, b) {
  if (is.list(a)) {
    a <- a[[1]]
  }
  # if (is.list(b)) {
  #   b <- b[[1]]
  # }
  cbind(intercept = 1, a) %**% b
}
