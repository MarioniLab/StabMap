#' imputeEmbedding
#'
#' imputeEmbedding
#'
#' @param assay_list assay_list
#' @param embedding embedding
#' @param reference reference
#' @param query query
#' @param neighbours neighbours
#' @param fun fun
#' @param ... ...
#'
#' @return matrix
#'
#' @examples
#'
#' @export
imputeEmbedding = function(assay_list,
                           embedding = NULL,
                           reference = Reduce(assay_list, colnames),
                           query = Reduce(assay_list, colnames),
                           neighbours = 5,
                           fun = mean,
                           ...) {

  # naive imputation given a (potentially batch corrected) StabMap embedding
  # input:
  # assay_list (typically the original input to StabMap)
  # StabMap embedding (if NULL then performs stabmap and passes along all the params)
  # which are the query cells
  # which are the reference cells
  # number of nearest neighbours
  # combining function: mean by default
  # default behaviour is to output a smoothed assay_list object

  require(BiocNeighbors)

  # given the embedding calculate the nearest neighbours

  has_reference = lapply(assay_list, function(x) any(reference %in% colnames(x)))

  imputed_list = list()

  for (assayName in names(assay_list)) {

    if (!has_reference[[assayName]]) next

    assayMat = assay_list[[assayName]]

    referenceCells = intersect(reference, colnames(assayMat))

    knn_out = queryNamedKNN(embedding[referenceCells, ], embedding[query,], neighbours)


    # cnames = colnames(assayMat)
    # require(Matrix)
    # imputedValues = apply(knn_out, 1, function(knnval)
    #     assayMat %*% (cnames %in% knnval)
    # )
    # rownames(imputedValues) <- rownames(assayMat)

    require(abind)
    imputedList = apply(knn_out, 2, function(knnval) {
      assayMat[,knnval]
    }, simplify = FALSE)
    imputedArray = abind(imputedList, along = 3)

    imputedMeans = apply(imputedArray, 1:2, fun)
    colnames(imputedMeans) <- rownames(knn_out)

    # ### example:
    # A <- array(c(rep(1,20), rep(2,20), rep(3,20)),dim = c(10,2,3))
    # B <- matrix(c(1:10), nrow = 2)
    # # multiply each A[,,i]%*%B
    #
    # C <- array(NA, dim=c(nrow(A), ncol(B), 3))
    # C[] <- apply(A, 3, function(x) x%*%B)
    # ###

    # want:
    # C = nrow(assayMat) x nrow(knn_out) x length(neighbours)
    # so therefore
    # B = assayMat
    # A is an array with dimensions
    # ncol(assayMat) x nrow(knn_out) x 5
    #
    # A_list = apply(knn_out, 2, function(knns) {
    #   cnames %in% knns
    # })

    # imputed_list[[assayName]] <- imputedValues / neighbours
    imputed_list[[assayName]] <- imputedMeans
  }

  return(imputed_list)

}
