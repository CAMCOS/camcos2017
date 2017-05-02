# Inout: data: sparse or dense form
#        n: number of dimensions to be (default is 10)
#        scale: the sclae of axis, in case the data is very sparse, try 1000, ..., 10^6, ...etc. if data is too sparse.
#


svdProcessing <- function(data,n=10, scale = 1){
  library(rgl) #for 3d plot
  library(irlba) # svd
  # SVD. set nv for different dimensions
  all.svd <- irlba(data, n)
  data.svd <- data %*% all.svd$v
  # 3d plot of data
  svdoutput <- plot3d(data.svd[,1], data.svd[,2], data.svd[,3], col="blue")
  snapshot3d("svdoutput") # sava as file svdout.pdf
  # 3d plot of data - rescale
  svdoutput2 <- plot3d(data.svd[,1], data.svd[,2], data.svd[,3], xlim = c(-10^(-scale), 10^(-scale)),
         ylim = c(-10^(-scale), 10^(-scale)),zlim = c(-10^(-scale), 10^(-scale)),col="blue")
  snapshot3d("svdoutput2")
  }