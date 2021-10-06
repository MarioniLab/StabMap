#' vectorSubset
#'
#' vectorSubset
#'
#' @param vec vec
#' @param mat mat
#'
#' @return matrix
#'
#' @keywords internal
vectorSubset = function(vec, mat) {
  # copied from SpatialUtils to avoid dependency
  # used for vectorised subsetting of a vector according to a matrix
  # within queryNamedKNN()

  # vec is a named vector
  # mat is a matrix containing the names or indices for which you want
  # to get the entries of vec

  vmat = c(mat)
  vvec = vec[vmat]

  vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
  colnames(vecmat) <- colnames(mat)
  rownames(vecmat) <- rownames(mat)

  return(vecmat)
}
