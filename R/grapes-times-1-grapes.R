#' Sorted matrix multiplication with intercept column
#'
#' This function first binds a column filled with 1s named \code{intercept} to
#' \code{a}, then performs rownames and colnames-aware (\code{\%**\%}) matrix
#' multiplication with \code{b}.
#'
#' @usage a \%*1\% b
#' @param a a matrix with rownames specified
#' @param b a matrix with colnames specified
#'
#' @return matrix
#'
#' @keywords internal
"%*1%" <- function(a, b) {
  if (is.list(a)) {
    a <- a[[1]]
  }
  cbind(intercept = 1, a) %**% b
}
