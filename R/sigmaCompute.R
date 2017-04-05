# D: distance
# k: k nearest neighbor, recommended value is log(n)

sigmaCompute <- function (D, k) 
{
    D = as.matrix(D)
    n = dim(D)[1]
    k = ifelse(k < 2, 2, k)
    D.sort = apply(D, 1, sort)
    dist.knn = D.sort[(k + 1), ]
    sigma = median(dist.knn)
    return(sigma)
}