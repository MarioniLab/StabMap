#' sparseCoef
#'
#' sparseCoef
#'
#' @param PCs_X PCs_X
#' @param PCs_X_genes PCs_X_genes
#'
#' @return matrix
#'
#' @examples
#'
#' @export
sparseCoef = function(PCs_X, PCs_X_genes) {
  require(glmnet)
  coefList_X = list()
  for (i in seq_len(ncol(PCs_X))) {
    # print(i)
    fit.cv = cv.glmnet(x = PCs_X_genes, y = PCs_X[,i], intercept = FALSE)
    fit = glmnet(x = PCs_X_genes, y = PCs_X[,i],
                 lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
    coefList_X[[i]] <- coef(fit)[-1]
  }
  coef_X = do.call(cbind,coefList_X)
  return(coef_X)
}
