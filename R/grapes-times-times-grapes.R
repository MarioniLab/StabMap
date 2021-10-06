#' Sorted matrix multiplication
#'
#' This function multiplies two matrices but first reorders the rows of the
#' second matrix to match the columns of the first matrix
#'
#' @usage X \%**\% Y
#' @param X a matrix with colnames specified.
#' @param Y a matrix with rownames specified. Alternatively, a list assumed to
#' contain two objects, a matrix with rownames specified, and a vector of
#' scaling values for subtraction.
#'
#' @return matrix
#'
#' @keywords internal
"%**%" <- function(X,Y) {
  # multiply two matrices but first reorder rows of the
  # second matrix

  # if Y is a list, then include a scaling subtraction
  # assumed given as second item

  if (is.list(Y)) {
    rmeans = Y[[2]]
    Y <- Y[[1]]
  } else {
    rmeans = rep(0, nrow(Y))
    names(rmeans) <- rownames(Y)
  }

  features = intersect(colnames(X), rownames(Y))
  XY = (X[,features] - rmeans[features]) %*% Y[features,]
  return(XY)
}
