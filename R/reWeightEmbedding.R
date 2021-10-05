#' Re-weight embedding
#'
#' Re-weights embedding according to given weights
#'
#' @param embedding Output from stabMap
#' @param weights weights vector
#' @param factor numeric value
#'
#' @return matrix
#'
#' @examples
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

  if (is.null(weights)) {
    weights = lapply(cols_split, function(x) 1)
  }

  norms = lapply(cols_split, function(cols) {
    sum(abs(embedding[,cols]))
  })

  norms_long = unsplit(norms, cols)
  weights_long = factor*unlist(weights)[cols]

  embedding_norm = t( (t(embedding) / norms_long) * weights_long)
  if (FALSE) {
    barplot(colSums(embedding_norm^2))
    boxplot(embedding_norm)
  }
  return(embedding_norm)
}
