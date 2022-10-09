
#' @importFrom matrixStats mean2 colMeans2
AddMean2draws <- function(x) {
  if (length(dim(x)) == 0) {
    c(x, mean2(x, refine = TRUE))
  } else {
    best <- colMeans2(x) |> array(dim = c(1, ncol(x)), dimnames = dimnames(x))
    rbind(x, best)
  }
}
