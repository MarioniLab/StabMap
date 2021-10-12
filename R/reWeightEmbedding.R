#' Re-weight StabMap embedding
#'
#' Re-weights embedding according to given weights for each reference dataset.
#' This gives more or less weighting to each contributing dataset and method
#' (PCA or LDA),
#'
#' @param embedding Joint embedding as output from stabMap.
#' @param weights (optional) named numeric vector giving relative weights for
#' each reference dataset.
#' @param factor numeric multiplicative value to offset near-zero values.
#'
#' @return matrix of same dimensions as `embedding`.
#'
#' @examples
#' set.seed(2021)
#' assay_list = mockMosaicData()
#' lapply(assay_list, dim)
#'
#' # specify which datasets to use as reference coordinates
#' reference_list = c("D1", "D3")
#'
#' # specify some sample labels to distinguish using linear discriminant
#' # analysis (LDA)
#' labels_list = list(
#' D1 = rep(letters[1:5], length.out = ncol(assay_list[["D1"]]))
#' )
#'
#' # stabMap
#' out = stabMap(assay_list,
#'               reference_list = reference_list,
#'               labels_list = labels_list,
#'               ncomponentsReference = 20,
#'               ncomponentsSubset = 20)
#'
#' # look at the scale of each component and discriminant
#' boxplot(out, las = 2, outline = FALSE)
#'
#' # re-weight embedding for less contribution from LDs and equal contribution
#' # from PCs of both references
#' out_reweighted = reWeightEmbedding(out, weights = c("D1_LD" = 0.5, "D1_PC" = 1, "D3_PC" = 1))
#'
#' # look at the new scale of each component and discriminant
#' boxplot(out_reweighted, las = 2, outline = FALSE)
#'
#' @export
reWeightEmbedding = function(embedding, weights = NULL, factor = 1e6) {

  # embedding is a cells x dimensions matrix
  # weights is an optional named list (names correspond to all before the underscore)
  # for weighting each embedding component
  # factor is a multiplicative value to avoid tiny numbers
  cols = as.character(interaction(gsub("_PC.*|_LD.*", "", colnames(embedding)),
                                  ifelse(grepl("_PC", colnames(embedding)), "PC", "LD"), sep = "_"))

  cols_split = split(colnames(embedding), cols)

  message(paste0(c("reweighting for references: ", names(cols_split)),
                 sep = " "))

  if (is.null(weights)) {
    weights = lapply(cols_split, function(x) 1)
  }

  norms = lapply(cols_split, function(cols) {
    sum(abs(embedding[,cols]))
  })

  norms_long = unsplit(norms, cols)
  weights_long = factor*unlist(weights)[cols]

  embedding_norm = t( (t(embedding) / norms_long) * weights_long)

  return(embedding_norm)
}
