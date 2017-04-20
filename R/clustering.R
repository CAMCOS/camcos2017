#' Spectral clustering
#'
#' @param Weights
#' @param method NJW, Ncut, DiffusionMap
#' @param k
#' @param t
#' @param sparse
#'
#' @return
#' @import RSpectra
#' @import Matrix
#' @export
#'
#' @examples
clustering <- function(Weights, method, k, t=NULL, sparse=T) {

    a <- Sys.time()

    require(RSpectra)
    require(Matrix)

    if (method == "NJW") {

        D <- rowSums(Weights)

        if( min(D) <= 0) {

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
        cluster.out <- kmeans(U, centers=k, nstart = 100)

    }


    else if (method == "Ncut") {

        n <- nrow(Weights)
        dvec_inv = 1/sqrt(rowSums(Weights))
        W_tilde = Matrix(rep(dvec_inv,n), ncol=n) * Weights * t(Matrix(rep(dvec_inv,n),ncol=n))
        W_tilde = (W_tilde+t(W_tilde))/2
        V <- eigs_sym(W_tilde, k, 'LM')$vector
        V = matrix(rep(dvec_inv,k-1), ncol = k-1)  * V[,2:k]

        cluster.out <- kmeans(V, k, nstart = 100)

    }


    else if (method == 'DiffusionMap'){

        if (is.null(t)) {

            stop("Specify a t value.") }

        n <- nrow(Weights)
        dvec_inv = 1/sqrt(rowSums(Weights))
        W_tilde = matrix(rep(dvec_inv,n), ncol=n) * Weights * t(matrix(rep(dvec_inv,n),ncol=n))
        W_tilde = (W_tilde+t(W_tilde))/2

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

    return(cluster.out$cluster)

}






