#' Spectral clustering
#'
#' @param Weights Matrix of similarity weights. Can be dense or sparse.
#' @param method NJW, Ncut, DiffusionMap.
#' @param k Number of clusters.
#' @param t Number of steps for the diffusion map.
#' @param sparse Logical. Indicates whether the input for the weights argument
#'   is sparse.
#' @param kmeans.method One of the following: "kmeans", "fuzzy", "poly.fuzzy"
#' @param m Fuzziness parameter. Defaults to m=2 for fuzzy; m=.5 for poly.fuzzy.
#'
#' @return Integer vector. Cluster assignments for each row of original data.
#' @import RSpectra
#' @import Matrix
#' @export
clustering <- function(Weights, method, k, t=NULL, sparse=TRUE,
                       kmeans.method="kmeans", m=NULL) {

  a <- Sys.time()

  #require(RSpectra)
  #require(Matrix)

  if (method == "NJW") {

    D <- rowSums(Weights)

    if( min(D) <= 0) { # if row has zero similarity, problems arise

      resp <- readline(prompt="One of your similarity rows has zero weight. Would you like to set
                       a 1 on the diagonal of the similarity? Type Y or N \n")
      if (resp == "Y" | resp == "y") {
        n <- which(D == 0)
        D[n] <- 1
      }

      else { stop("One of your rows has zero weight.") }
    }

    D <- Diagonal(n=nrow(Weights),(D^-.5))
    Z <- D %*% Weights %*% D

    if (sparse==F) {
      EZ <- eigs_sym(Z, k, 'LM')$vector }
    else {
      require(irlba)
      EZ <- partial_eigen(x=Z, n = k, symmetric = TRUE)$vectors
    }
    U <- EZ/sqrt(rowSums(EZ^2))

    # run kmeans in eigenspace:

    if (kmeans.method=="kmeans") {
      cluster.out <- kmeans(U, centers=k, nstart = 100)
    }

    else if (kmeans.method=="fuzzy") {
      require(fclust)
      if (is.null(m)) {
        m <- 2
      }
      cluster.out <- FKM(X=U,k=k, m=m, RS=10)
    }

    else if (kmeans.method=="poly.fuzzy") {
      require(fclust)
      if (is.null(m)) {
        m <- .5
      }
      cluster.out <- FKM.pf(X=U,k=k, b=m, RS=10)
    }
}


  else if (method == "Ncut") {

    n <- nrow(Weights)
    dvec_inv = 1/sqrt(rowSums(Weights))
    W_tilde = Matrix(rep(dvec_inv,n), ncol=n) * Weights * t(Matrix(rep(dvec_inv,n),ncol=n))
    W_tilde = (W_tilde+t(W_tilde))/2

    # diag(dvec_inv) %*% Weights %*% diag(dvec_inv) ?
    # why the average part?

    V <- eigs_sym(W_tilde, k, 'LM')$vector
    V = matrix(rep(dvec_inv,k-1), ncol = k-1) * V[,2:k]
    V = V / (matrix(rep(sqrt(rowSums(V^2)),k-1),ncol=k-1))
    cluster.out <- kmeans(V, k, nstart = 100)

  }


  else if (method == 'DiffusionMap'){

    if (is.null(t)) {

      stop("Specify a t value.") }

    n <- nrow(Weights)
    dvec_inv = 1/sqrt(rowSums(Weights))
    W_tilde = matrix(rep(dvec_inv,n), ncol=n) * Weights * t(matrix(rep(dvec_inv,n),ncol=n))
    W_tilde = (W_tilde+t(W_tilde))/2

    # diag(dvec_inv) %*% Weights %*% diag(dvec_inv) ?
    # why the average part?

    V <- eigs_sym(W_tilde, k, 'LM')$vector
    lambda <- eigs_sym(W_tilde, k, 'LM')$value
    V_inv = 1/sqrt((rowSums(V[,2:k]^2)))
    V <- matrix(rep(V_inv,k-1), ncol=k-1) * V[,2:k]
    V = matrix(rep(dvec_inv,k-1), ncol = k-1)  * V
    V = (matrix(rep(lambda[2:k], each=n), ncol=k-1)^t )* V

    cluster.out <- kmeans(V, k, nstart = 100)

  }

  b <- Sys.time()
  time.elapsed <- b - a
  print(time.elapsed)

  return(cluster.out)

}






