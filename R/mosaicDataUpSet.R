#' mosaicDataUpSet
#'
#' Plots feature overlaps of mosaic data as an UpSet plot.
#'
#' @param assay_list a list of data matrices with rownames (features) specified.
#' @param plot logical (default FALSE) whether the UpSet plot should be printed.
#' @param ... further arguments passed to `upset` from the `UpSetR` package.
#'
#' @return UpSet object displaying degree of overlap of rownames (features)
#' among each of the data matrices in \code{assay_list}. Set bars correspond to
#' the number of cells/samples present in each data matrix.
#'
#' @examples
#' set.seed(2021)
#' assay_list = mockMosaicData()
#' lapply(assay_list, dim)
#' mosaicDataUpSet(assay_list)
#'
#' # additional arguments from UpSetR::upset()
#' mosaicDataUpSet(assay_list, empty.intersections = TRUE)
#'
#' @export
mosaicDataUpSet = function(assay_list, plot = FALSE, ...) {
  require(UpSetR)

  df = as.data.frame(1 * do.call(
    cbind,
    lapply(
      assay_list,
      function(x) Reduce(union, lapply(assay_list, rownames)) %in%
        rownames(x))))

  df_cells =  as.data.frame(1 * do.call(
    cbind, lapply(
      assay_list,
      function(x) Reduce(union, lapply(assay_list, colnames)) %in%
        colnames(x))))

  g <- upset(df,...)

  g0 <- suppressMessages(upset(df_cells,...))
  g$Sizes <- g0$Sizes

  if (plot) {
    print(g)
  }
  return(g)
}
