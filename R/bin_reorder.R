#' "%**%"
#'
#' "%**%"
#'
#' @param X matrix
#' @param Y matrix
#'
#' @return matrix
#'
#' @examples
#'
#' @export
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
