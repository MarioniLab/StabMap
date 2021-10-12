#' mockMosaicData
#'
#' Mock up a mosaic data list using simulated data, for use in documentation
#' examples.
#'
#' @param names character vector of mock datasets.
#' @param ncells integer vector of cells in each mock dataset.
#' @param ngenes list containing integer vectors of features measured in each
#' mock dataset.
#' @param fun name of function to simulate data, default "rnorm".
#' @param ... further arguments passed to `fun`.
#'
#' @return assay_list a list of data matrices with rownames (features)
#' specified.
#'
#' @examples
#' set.seed(2021)
#' assay_list = mockMosaicData()
#' lapply(assay_list, dim)
#'
#' # simulate data from another distribution
#' assay_list = mockMosaicData(fun = "rnbinom", size = 5, prob = 0.5)
#' lapply(assay_list, dim)
#'
#' @export
mockMosaicData = function(names = c("D1", "D2", "D3"),
                          ncells = c(50, 50, 50),
                          ngenes = list(1:150,76:225,151:300),
                          fun = "rnorm",
                          ...) {

  assay_list = mapply(function(name, ncell, ngene)
    matrix(get(fun)(ncell*length(ngene), ...),
           nrow = length(ngene),
           ncol = ncell,
           dimnames = list(
             paste0("gene_", ngene),
             paste0(name, "_cell_", seq_len(ncell))
           )),
    names, ncells, ngenes, SIMPLIFY = FALSE)

  return(assay_list)
}
