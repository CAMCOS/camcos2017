

#' Calculate similarity matrices
#'
#' Note: this function sets the self-similarity diagonal = 0. For certain
#' applications, diagonal = 1 may be more appropriate.
#'
#' @param data Matrix.
#' @param method One of the following: "correlation", "cosine", "dotproduct", "Gaussian".
#' @param rowscaling One of the following: "L2", "L1".
#' @param colscaling Single option: "standardize" (center and divide by variance).
#' @param sigma Variance for Gaussian Kernel method.
#' @param centers Proportion of columns to be sampled, if desired (e.g. centers = .25).
#' @param seed For sampling reproducibility.
#' @param distance For use with Gaussian Kernel; options = "JSdivergence", "L1", "L2"
#' @param sparse Logical. Is the data sparse or not.
#'
#' @return Similarity matrix.
#' @import Matrix
#' @export
similarity <- function(data, method, rowscaling = NULL, colscaling = NULL,
                       sigma = NULL, centers = NULL, seed = NULL,
                       distance = NULL, sparse = TRUE) {

  a <- Sys.time()


  # data <- Matrix(data)

  ##### Column scaling #####

  if (!(is.null(colscaling))) {
    if (colscaling == "standardize") {
      data <- apply(data, 2, scale) }

    else {stop("Pick a valid column scaling.")}

  }

  #####
  #####
  ##### Row scaling #####

  if (!(is.null(rowscaling))) {

    if (rowscaling == "L2") {
      data <- data/sqrt(rowSums(data^2)) }

    else if (rowscaling == "L1") {
      data <- data/rowSums(data) }

    else {stop("Pick a valid row scaling.")}

  }

  #####
  #####
  ##### Distance -> Gaussian similarity #####

  if (method == "Gaussian") {

    # Calculate a distance metric.

    if (distance == "JSdivergence") {

      jsdiv <- function(P){
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
        return(D)
      }

      dist <- jsdiv(data) }

    else if (distance == "L2") {

      dist <- as.matrix(dist(data)) }

    else if (distance == "L1") {

      dist <- as.matrix(dist(data, method = "manhattan")) }

    else {stop("Pick a valid distance.")}


    # Convert distance to similarity

    if (is.null(sigma)) {stop("Choose a sigma value.")}

    else { Similarity <- exp(-1 * dist^2 / (2*sigma)) }

  }

  #####
  #####
  ##### Full data similarity matrix #####

  # dense matrix:

  if (sparse == F) {

    if (is.null(centers)) {

      if (method == "correlation") {

        centeredcolumns <- t(t(data)-colMeans(data)) # center the data by column

        rowvar <- rowSums(centeredcolumns^2) # store the variances for each row (where colMeans=0)

        Cov.matrix <- tcrossprod(centeredcolumns) # calculate the covariance matrix (dot product all rows)

        # corr = cov(x,y) / sqrt(var(x)var(y)
        Corr.true <- Cov.matrix/sqrt(rowvar)
        Similarity <- t(t(Corr.true)/sqrt(rowvar))

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(0, nrow(Similarity))

      }

      else if (method == "corr.hack") {

        # center the data by column
        centeredcolumns <- t(t(data)-colMeans(data))
        # # store the variances for each row (where colMeans=0)
        # rowvar <- rowSums(centeredcolumns^2)
        # calculate the covariance matrix (dot product all rows)
        Cov.matrix <- tcrossprod(centeredcolumns)
        # calculate the length of each row (eventually scale by row&col)
        cov.rowsums <- rowSums(Cov.matrix^2)
        Corr.hack <- Cov.matrix/sqrt(cov.rowsums)
        Similarity <- t(t(Corr.hack)/sqrt(cov.rowsums))
        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(0, nrow(Similarity))

      }

      else if (method == "cosine") {

        rowlength <- rowSums(data^2)
        Dot.prods <- tcrossprod(data)
        Cosines <- Dot.prods/sqrt(rowlength)
        Similarity <- t(t(Cosines)/sqrt(rowlength))
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(0, nrow(Similarity))

      }

      else if (method == "dotproduct") {

        Similarity <- tcrossprod(data)
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(0, nrow(Similarity))

      }

    }

    #####
    #####
    ##### Random centers similarity matrix #####

    else {

      if (is.numeric(seed)) { set.seed(seed) }

      if (method == "correlation") {

        centeredcolumns <- data-colMeans(data) # center the data by column

        rowvar.full <- rowSums(centeredcolumns^2) # store the variances for each row (where colMeans=0)

        # rcenters is rxN matrix, r = kcenters % of rows, sampled randomly
        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(centeredcolumns[sample(nrow(data),kcenters,replace = FALSE), ])

        rowvar.centers <- rowSums(rcenters^2)

        dotprod.centers <- tcrossprod(centeredcolumns,rcenters)

        # corr = cov(x,y) / sqrt(var(x)var(y)
        Corr.true <- dotprod.centers/sqrt(rowvar.full)
        Corr.centers <- t(t(Corr.true)/sqrt(rowvar.centers))

        Similarity <- tcrossprod(Corr.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

      else if (method == "corr.hack") {

        centeredcolumns <- data-colMeans(data) # center the data by column

        # rcenters is rxN matrix, r = kcenters % of rows, sampled randomly
        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(centeredcolumns[sample(nrow(data),kcenters,replace = FALSE), ])

        dotprod.centers <- tcrossprod(centeredcolumns,rcenters)

        rowlengths <- rowSums(dotprod.centers)
        collengths <- colSums(dotprod.centers)

        # hack = cov(x,y) / sqrt(length(x)length(y))
        Corr.hack <- dotprod.centers/sqrt(rowlengths)
        Corrhack.centers <- t(t(Corr.hack)/sqrt(collengths))

        Similarity <- tcrossprod(Corrhack.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

      else if (method == "cosine") {

        rowlengths.full <- rowSums(data^2) # store the variances for each row (where colMeans=0)

        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(data[sample(nrow(data),kcenters,replace = FALSE), ])

        rowlengths.centers <- rowSums(rcenters^2)

        dotprod.centers <- tcrossprod(data,rcenters)

        # cos = <x,y> / sqrt(length(x)length(y))
        Cosine <- dotprod.centers/sqrt(rowlengths.full)
        Cosine.centers <- t(t(Cosine)/sqrt(rowlengths.centers))

        Similarity <- tcrossprod(Cosine.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

      else if (method == "dotproduct") {

        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(data[sample(nrow(data),kcenters,replace = FALSE), ])

        dotprod.centers <- tcrossprod(data,rcenters)

        Similarity <- tcrossprod(dotprod.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

    }

    #####
    #####

  }


  # sparse matrix:

  else {

    if (is.null(centers)) {

      if (method == "correlation") {

        # centeredcolumns <- data
        #
        # rowvar <- rowSums(centeredcolumns^2) # store the variances for each row (where colMeans=0)
        #
        # Cov.matrix <- tcrossprod(centeredcolumns) # calculate the covariance matrix (dot product all rows)
        #
        # # corr = cov(x,y) / sqrt(var(x)var(y)
        # Corr.true <- Cov.matrix/sqrt(rowvar)
        # Similarity <- t(t(Corr.true)/sqrt(rowvar))
        #
        # # scale the matrix to (0,1), excluding the diagonal
        # diag(Similarity) <- rep(0,nrow(Similarity))
        # Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        # diag(Similarity) <- rep(0, nrow(Similarity))

        stop("Correlation does not work on a sparse matrix. Try using a dense matrix instead.")

      }

      else if (method == "corr.hack") {

        # centeredcolumns <- data
        # # # store the variances for each row (where colMeans=0)
        # # rowvar <- rowSums(centeredcolumns^2)
        # # calculate the covariance matrix (dot product all rows)
        # Cov.matrix <- tcrossprod(centeredcolumns)
        # # calculate the length of each row (eventually scale by row&col)
        # cov.rowsums <- rowSums(Cov.matrix^2)
        # Corr.hack <- Cov.matrix/sqrt(cov.rowsums)
        # Similarity <- t(t(Corr.hack)/sqrt(cov.rowsums))
        # # scale the matrix to (0,1), excluding the diagonal
        # diag(Similarity) <- rep(0,nrow(Similarity))
        # Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        # diag(Similarity) <- rep(0, nrow(Similarity))

        stop("Correlation does not work on a sparse matrix. Try using a dense matrix instead.")

      }

      else if (method == "cosine") {

        rowlength <- rowSums(data^2)
        data <- data/sqrt(rowlength)    # cosine normalizes each vector to unit length
        Similarity <- tcrossprod(data)
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(0, nrow(Similarity))

      }

      else if (method == "dotproduct") {

        Similarity <- tcrossprod(data)
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(0, nrow(Similarity))

      }

    }

    #####
    #####
    ##### Random centers similarity matrix #####

    else {

      if (is.numeric(seed)) { set.seed(seed) }

      if (method == "correlation") {

        centeredcolumns <- data-colMeans(data) # center the data by column

        rowvar.full <- rowSums(centeredcolumns^2) # store the variances for each row (where colMeans=0)

        # rcenters is rxN matrix, r = kcenters % of rows, sampled randomly
        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(centeredcolumns[sample(nrow(data),kcenters,replace = FALSE), ])

        rowvar.centers <- rowSums(rcenters^2)

        dotprod.centers <- tcrossprod(centeredcolumns,rcenters)

        # corr = cov(x,y) / sqrt(var(x)var(y)
        Corr.true <- dotprod.centers/sqrt(rowvar.full)
        Corr.centers <- t(t(Corr.true)/sqrt(rowvar.centers))

        Similarity <- tcrossprod(Corr.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

      else if (method == "corr.hack") {

        centeredcolumns <- data-colMeans(data) # center the data by column

        # rcenters is rxN matrix, r = kcenters % of rows, sampled randomly
        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(centeredcolumns[sample(nrow(data),kcenters,replace = FALSE), ])

        dotprod.centers <- tcrossprod(centeredcolumns,rcenters)

        rowlengths <- rowSums(dotprod.centers)
        collengths <- colSums(dotprod.centers)

        # hack = cov(x,y) / sqrt(length(x)length(y))
        Corr.hack <- dotprod.centers/sqrt(rowlengths)
        Corrhack.centers <- t(t(Corr.hack)/sqrt(collengths))

        Similarity <- tcrossprod(Corrhack.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

      else if (method == "cosine") {

        rowlengths.full <- rowSums(data^2) # store the variances for each row (where colMeans=0)

        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(data[sample(nrow(data),kcenters,replace = FALSE), ])

        rowlengths.centers <- rowSums(rcenters^2)

        dotprod.centers <- tcrossprod(data,rcenters)

        # cos = <x,y> / sqrt(length(x)length(y))
        Cosine <- dotprod.centers/sqrt(rowlengths.full)
        Cosine.centers <- t(t(Cosine)/sqrt(rowlengths.centers))

        Similarity <- tcrossprod(Cosine.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

      else if (method == "dotproduct") {

        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(data[sample(nrow(data),kcenters,replace = FALSE), ])

        dotprod.centers <- tcrossprod(data,rcenters)

        Similarity <- tcrossprod(dotprod.centers)

        # scale the matrix to (0,1), excluding the diagonal
        diag(Similarity) <- rep(0,nrow(Similarity))
        Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
        diag(Similarity) <- rep(1, nrow(Similarity))

      }

    }

    #####
    #####

  }  #

  b <- Sys.time()
  time.elapsed <- b - a

  print(time.elapsed)

  # Output
  Similarity

}
