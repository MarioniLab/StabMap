#' plotFeatureOverlaps
#'
#' plots feature overlaps
#'
#' @param assay_list assay list
#'
#' @return upset plot
#'
#' @examples
#'
#' @export
plotFeatureOverlaps = function(assay_list) {
  require(UpSetR)
  g = upset(as.data.frame(1*do.call(cbind, lapply(assay_list, function(x) Reduce(union, lapply(assay_list, rownames)) %in% rownames(x)))))
  print(g)
  return(g)
}
