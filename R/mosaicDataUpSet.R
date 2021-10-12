#' mosaicDataUpSet
#'
#' Plots feature overlaps of mosaic data as an UpSet plot.
#'
#' @param assay_list a list of data matrices with rownames (features) specified.
#' @param plot logical (default TRUE) whether the UpSet plot should be printed.
#' @param ... further arguments passed to `upset` from the `UpSetR` package.
#'
#' @return UpSet object displaying degree of overlap of rownames (features)
#' among each of the data matrices in \code{assay_list}.
#'
#' @examples
#' assay_list = mockMosaicData()
#' mosaicDataUpSet(assay_list)
#'
#' # additional arguments from UpSetR::upset()
#' mosaicDataUpSet(assay_list, empty.intersections = TRUE)
#'
#' @export
mosaicDataUpSet = function(assay_list, plot = TRUE, ...) {
  require(UpSetR)
  g = upset(as.data.frame(1*do.call(cbind, lapply(assay_list, function(x)
    Reduce(union, lapply(assay_list, rownames)) %in% rownames(x)))),
    ...)
  if (plot) {
    print(g)
  }
  return(g)
}
