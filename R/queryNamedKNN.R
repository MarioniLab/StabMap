#' queryNamedKNN
#'
#' queryNamedKNN
#'
#' @param coords_reference coords_reference
#' @param coords_query coords_query
#' @param k k
#'
#' @return matrix
#'
#' @keywords internal
queryNamedKNN = function(coords_reference, coords_query, k) {
  # used in imputeEmbedding()

  require(BiocNeighbors)

  knn = queryKNN(
    coords_reference,
    coords_query,
    k = k, get.distance = FALSE)$index
  rownames(knn) <- rownames(coords_query)
  knn_name = vectorSubset(rownames(coords_reference), knn)

  return(knn_name)
}
