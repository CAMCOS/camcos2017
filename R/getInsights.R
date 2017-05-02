#' Determine the Top-N Keywords Per Cluster
#'
#' @param cluster Data from one of the clusters.
#' @param vocab Column names, like vocabulary, websites. Vocab should be one to
#'   one with column indices.
#' @param n Top n features you would like to see.
#' @param plot Logical. To produce plot in addition to outputting top-n keywords and
#'   maginitudes.
#' @param file Filepath and PNG file name. Note: plot must be set to TRUE to save a
#'   file.
#'
#' @return Barchart of the top-n influential keywords in the cluster.
#' @import RSpectra
#' @import lattice
#' @export
getInsights <- function(cluster, vocab, n, plot = FALSE, file = NULL){
    #library(Rspectra)
    svd.out <- svds(cluster, 4)
    v <- svd.out$v #dim(v) 61066     4

    b_clr <- c("steelblue", "darkred")
    key <- simpleKey(rectangles = TRUE, space = "top", points=FALSE,
      text=c("Positive", "Negative"))
    key$rectangles$col <- b_clr

    v1_top <- order(abs(v[,1]), decreasing = T)[1:n]
    v1_top_po <- v1_top[which(v[,1][v1_top] > 0)] # positive values
    v1_top_ne <- v1_top[which(v[,1][v1_top] < 0)] # negative values
    v1_top_t <- c(v1_top_po,v1_top_ne)
    v1_topn <- v[,1][v1_top_t]
    v1_top.words <- as.matrix(vocab[v1_top_t])

    if (plot) {
      b1 <- barchart(as.table(v1_topn),
                     main="First column",
                     horizontal=FALSE, col=ifelse(v1_topn > 0,
                                                  b_clr[1], b_clr[2]),
                     ylab="Impact value",
                     scales=list(x=list(rot=55, labels=v1_top.words, cex=0.9)),
                     key = key)
      if (!is.null(file)) {
        png(file)
        print(b1)
        dev.off()
      } else {
        print(b1)
      }
    }

    list(magnitude = v1_topn, keywords = as.character(v1_top.words))

}

