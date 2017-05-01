# Input   cluster: data from one of the clusters
#          vocab: column names, like vocabulary, websites...
#		   n : top n features you would like to see
# Note: vocab should be one to one with column indices.

getInsights <- function(cluster, vocab, n){
  library(Rspectra)
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

b1 <- barchart(as.table(v1_topn),
  main="First column",
  horizontal=FALSE, col=ifelse(v1_topn > 0, 
      b_clr[1], b_clr[2]),
  ylab="Impact value", 
  scales=list(x=list(rot=55, labels=v1_top.words, cex=0.9)),
  key = key)
b1
}

