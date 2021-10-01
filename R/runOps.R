#' runOps
#'
#' runOps
#'
#' @param obj obj
#' @param ops ops
#' @param leftToRight leftToRight
#'
#' @return matrix
#'
#' @examples
#'
#' @export
runOps = function(obj, ops, leftToRight = TRUE) {
  if (leftToRight) {
    out = obj[[1]]
    for (i in 1:length(ops)) {
      out <- get(ops[[i]])(out, obj[[i+1]])
    }
  } else {
    out = obj[[length(obj)]]
    for (i in length(ops):1) {
      out <- get(ops[[i]])(obj[[i]], out)
    }
  }
  return(out)
}
