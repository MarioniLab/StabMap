#' "%projpred%"
#'
#' "%projpred%"
#'
#' @param a a
#' @param b b
#'
#' @return matrix
#'
#' @examples
#'
#' @export
"%projpred%" <- function(a, b) {
  # a is a matrix
  # b is a list of a 1. projection matrix
  # and 2. a model to pass through predict
  # if b is not a list, then just do normal
  # matrix multiplication
  # alternatively, if b is already an lda
  # object then just perform the prediction

  if (class(b)[1] == "lda") {
    features = rownames(b$scaling)
    am = MASS:::predict.lda(b, newdata = a[,features])$x
    return(am)
  }

  if (!is.list(b)) {
    return(a %*% b)
  }

  ab = a %*% b[[1]]
  if (class(b[[2]]) == "lda") {
    am = MASS:::predict.lda(b[[2]], newdata = a)$x
  }
  if (class(b[[2]]) == "svm") {
    am = attr(predict(b[[2]], newdata = a, decision.values = TRUE), "decision.values")
  }

  return(cbind(ab,am))

}
