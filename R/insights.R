#' Determine the Top-N Keywords Per Cluster
#'
#' @param data
#' @param labels Integer vector. Labels output by `clustering()`
#' @param cluster Integer. Cluster label for which to calculate top keywords.
#' @param top_n Integer. Specificies the number of top keywords to show.
#' @param plot Logical. If TRUE, a barplot of the keyword importance is returned.
#'
#' @return List of 2. First element: keywords in decreasing order of importance
#'   Second element: maginitude of importance associated with the keyword.
#' @export
insights <- function(data, labels, cluster, top_n = 10, plot = TRUE, file = NULL) {

}
