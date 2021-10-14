#' Project and/or predict data using feature weights or a LDA model object
#'
#' This function takes a data matrix a and, depending on the class of b,
#' projects the data using feature weights, or predicts new values using
#' linear discriminant analysis (LDA) model object, or both.
#'
#' @usage a \%projpred\% b
#' @param a a matrix with colnames specified
#' @param b a matrix with rownames specified, or a lda model object, or a
#' list containing a matrix and/or a lda model object.
#'
#' @return matrix
#'
#' @keywords internal
"%projpred%" <- function(a, b) {
  # a is a matrix
  # b is a list of a 1. projection matrix
  # and 2. a model to pass through predict
  # if b is not a list, then just do normal
  # matrix multiplication
  # alternatively, if b is already an lda
  # object then just perform the prediction

  if (is(b, "lda")) {
    require(MASS)
    features = rownames(b$scaling)
    am = MASS:::predict.lda(b, newdata = a[,features])$x
    return(am)
  }

  if (!is.list(b)) {
    return(a %*% b)
  }

  ab = a %*% b[[1]]
  if (is(b[[2]], "lda")) {
    require(MASS)
    am = MASS:::predict.lda(b[[2]], newdata = a)$x
  }
  if (is(b[[2]], "svm")) {
    am = attr(predict(b[[2]], newdata = a, decision.values = TRUE), "decision.values")
  }

  return(cbind(ab,am))

}
