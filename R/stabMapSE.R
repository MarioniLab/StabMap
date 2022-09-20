#' Stabilised mosaic single cell data integration using unshared features
#'
#' stabMapSE performs StabMap with SummarizedExperiment type
#' (SingleCellExperiment, SpatialExperiment) objects as input.
#'
#' @param ... Set of SummarizedExperiment, SingleCellExperiment,
#' SpatialExperiment, etc that should be named.
#' @param assays Named character vector of assays to be extracted from the
#' objects.
#' @param args List of arguments to be passed on to `stabMap`
#'
#' @return matrix containing common embedding with rows corresponding to cells,
#' and columns corresponding to PCs or LDs for reference dataset(s).
#'
#' @examples
#' set.seed(2021)
#' library(SingleCellExperiment)
#'
#' sce_1 = SingleCellExperiment(assays = list(logcounts = mockMosaicData(names = "D1", ncells = 1000, ngenes = list(1:500))[[1]]))
#' sce_2 = SingleCellExperiment(assays = list(logcounts = mockMosaicData(names = "D2", ncells = 1000, ngenes = list(251:750))[[1]]))
#' sce_3 = SingleCellExperiment(assays = list(counts = mockMosaicData(names = "D3", ncells = 1000, ngenes = list(500:750))[[1]]))
#'
#' out = stabMapSE(D1 = sce_1, D2 = sce_2)
#'
#' # pass on additional parameters to stabMap()
#' # e.g. change number of components of reference
#' out = stabMapSE(D1 = sce_1, D2 = sce_2, args = list(ncomponentsReference = 20))
#' dim(out)
#'
#' # pull out different assay names:
#' out = stabMapSE(D2 = sce_2, D3 = sce_3, assays = c("D2" = "logcounts", "D3" = "counts"))
#'
#' @export
stabMapSE = function(...,
                     assays = "logcounts",
                     args = list()
                     ) {

  SE_list = list(...)

  if (identical(assays,"logcounts")) {
    assays = setNames(rep("logcounts", length(SE_list)), names(SE_list))
  }

  assay_list = mapply(assay, SE_list, assays[names(SE_list)], SIMPLIFY = FALSE)

  all_embeddings = do.call(stabMap, c(list(assay_list = assay_list), args))

  return(all_embeddings)
}
