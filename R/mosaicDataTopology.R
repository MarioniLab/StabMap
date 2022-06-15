#' mosaicDataTopology
#'
#' Generate mosaic data topology network as an igraph object.
#'
#' @param assay_list a list of data matrices with rownames (features) specified.
#'
#' @return igraph weighted network with nodes corresponding to
#' \code{assay_list} elements, and edges present if the matrices share at
#' least one rowname. Edge weights correspond to the number of shared
#' rownames among data matrices.
#'
#' @examples
#' set.seed(2021)
#' assay_list = mockMosaicData()
#' mdt = mosaicDataTopology(assay_list)
#' mdt
#' plot(mdt)
#'
#' @export
mosaicDataTopology = function(assay_list) {
  require(igraph)
  # given a list of assays, generate a
  # network relating the datasets to each other
  # in terms of number of shared features
  # (rownames)

  datasets = names(assay_list)

  pairs = t(combn(datasets, 2))

  edge_weights = apply(pairs, 1, function(x) {
    length(Reduce(intersect, lapply(assay_list[x], rownames)))
  })

  pairs_overlapping = pairs[edge_weights != 0,, drop = FALSE]
  edge_weights_overlapping = edge_weights[edge_weights != 0]

  g = graph.edgelist(pairs_overlapping, directed = FALSE)
  E(g)$weight <- edge_weights_overlapping

  g = graph.edgelist(pairs_overlapping, directed = FALSE)
  sd = setdiff(datasets, V(g)$name)
  if (length(sd) > 0) {
    g <- add_vertices(g, length(sd), name = sd)
  }

  if (components(g)$no != 1) {
    message("feature network is not connected, features must overlap via rownames for StabMap to run")
  }

  # add some aesthetic attributes to the network
  V(g)$frame.color = "white"
  V(g)$color = "white"
  V(g)$label.color = "black"
  V(g)$label.family = "sans"

  return(g)
}
