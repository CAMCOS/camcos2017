#Because the KL Divergence is not symmetric, and gives a value of infinity for 
#distributions that lack the same support values, the Jensen-Shannon Divergence can be better
#what this does, is compute a M Matrix that is the average of P and Q, then perform KL(P,M) and KL (Q,M)
#then averages those values.
jsdiv <- function(P){
  if(!is.matrix(P)) stop("Input must be a matrix of row Probability Distributions")
  nrows <- length(P[,1])
  ncols <- length(P[1,])
  D <- matrix(rep.int(0, nrows ** 2), nrow = nrows)
  P[is.nan(P)] <- 0
  for(i in 2:nrows){
    p.row <- P[i,]
    for(j in 1:i-1){
      q.row <- P[j,]
      m.row <- 1/2 * (p.row + q.row)
      D[i,j] <- D[j,i] <- (1/2 * sum(p.row * log(p.row/m.row), na.rm =TRUE) + 1/2 * sum(q.row * log(q.row/m.row), na.rm = TRUE))
    }
  }
  D
}