#' Run a sequence of binary operations
#'
#' @param obj list of objects.
#' @param ops list of operations (length should be 1 less than `obj`).
#' @param leftToRight logical whether operations should be performed in order
#' from left to right (default), or right to left.
#'
#' @return matrix or array output of the sequence of binary operations
#'
#' @keywords internal
.runOps = function(obj, ops, leftToRight = TRUE) {
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
