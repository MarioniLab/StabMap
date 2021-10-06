#' mosaicDataUpSet
#'
#' plots feature overlaps as an UpSet plot
#'
#' @param assay_list a list of matrices with rownames (features) specified
#'
#' @return UpSet object displaying degree of overlap of rownames (features)
#' among each of the data matrices in \code{assay_list}.
#'
#' @examples
#'
#' @export
mosaicDataUpSet = function(assay_list, plot = TRUE) {
  # previously named plotFeatureOverlaps
  require(UpSetR)
  g = upset(as.data.frame(1*do.call(cbind, lapply(assay_list, function(x)
    Reduce(union, lapply(assay_list, rownames)) %in% rownames(x)))))
  if (plot) {
    print(g)
  }
  return(g)
}
