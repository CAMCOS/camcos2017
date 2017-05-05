#' Dimensionality Reduction Via SVD and 3-D Plotting
#'
#' @param data Matrix object. Accepts sparse or dense form.
#' @param n Number of dimensions to be reduced to (default is 10)
#' @param scale The sclae of axis, in case the data is very sparse. If data is
#'   too sparse, try 1000, ..., 10^6, etc.
#' @param plot Logical. Return 3-D plots
#' @param file Full path to filename without file type, which will be PNG.
#'
#' @return Matrix object with the same number of rows as the original input
#'   data, and a reduced number of columns equal to the input argument \code{n}.
#' @import rgl
#' @import irlba
#' @export
svdProcessing <- function(data, n=10, scale = 1, plot = FALSE, file = NULL) {

  # SVD. set nv for different dimensions
  all.svd <- irlba(data, n)
  data.svd <- data %*% all.svd$v

  if (plot) {
    svdoutput <- plot3d(data.svd[,1], data.svd[,2], data.svd[,3], col="blue")
    if (!is.null(file)) {
      snapshot3d(file)
    }
    svdoutput2 <- plot3d(data.svd[,1], data.svd[,2], data.svd[,3], xlim = c(-10^(-scale), 10^(-scale)),
                         ylim = c(-10^(-scale), 10^(-scale)),zlim = c(-10^(-scale), 10^(-scale)),col="blue")
    if (!is.null(file)) {
      snapshot3d(paste(file, "_scaled.png", sep = ""))
    }
  }

  data.svd
}
