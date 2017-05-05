#' Column weighting
#'
#' @param data Matrix of counts.
#' @param weightfunction One of the following: beta", "step", "linear", "IDF",
#'   "IDF^2".
#' @param sparseinput Logical. Is your data in sparse (long) format? (Long
#'   format = Nx3 matrix-like object of nonzero entries)
#' @param par1 Beta parameters or step function boundaries.
#' @param par2 Beta parameters or step function boundaries.
#' @param mode Linear function max (desired density max weight, e.g. 1/k).
#' @param binary Logical. Convert data to binary?
#' @param convertsparse Logical. If your matrix is dense, convert to sparse
#' @param lower Integer. Lower bound for column sum in order to retain column;
#'   columns with sums below this threshold will be removed.
#' @param upper Integer. Upper bound for column sum in order to retain column;
#'   columns with sums above this threshold will be removed.
#'
#' @return Processed version of original input matrix.
#' @seealso \code{\link{similarity}}
#' @import Matrix
#' @export
colweights <- function (data, weightfunction, sparseinput,
                        par1=NULL, par2=NULL, mode=NULL,
                        binary=TRUE,
                        convertsparse=TRUE,
                        lower=2, upper=NULL) {

  # require(Matrix)

  #####
  ##### construct Matrix object for use with Matrix package #####

  if (sparseinput==T) {    # given a sparse matrix - convert to Matrix class

    if (is.matrix(data)) {
      if (binary==T) {
        data <- sparseMatrix(i=data[,1], j=data[,2], x=rep(1, nrow(data)))
      }

      else {
        data <- sparseMatrix(i=data[,1], j=data[,2], x=data[,3])
      }
    }

    else if (is.data.frame(data)) {

      if (binary==T) {
        data <- sparseMatrix(i=data[,1], j=data[,2], x=rep(1, nrow(data)))
      }

      else {
        data <- sparseMatrix(i=data[,1], j=data[,2], x=data[,3])
      }
    }

    else if (binary==T) {

      data[data>0]<-1

    }

  }

  else{   # given a dense matrix (sparseinput == F)

    if(binary==TRUE) {

      data[data>0]<-1

    }

    if(convertsparse==TRUE) {   # convert dense matrix to sparse matrix
      if (is.matrix(data)) {
        data <- Matrix(data, sparse=T)
      }

      else if (is.data.frame(data)) {
        data <- as.matrix(data)
        data <- Matrix(data)
      }
    }

    else {  # keep data in dense format, but convert to class = Matrix
      data <- as.matrix(data)
      data <- Matrix(data)
    }

  }



  #####
  #####
  ##### Find density proportion of each column #####

  weightfunction <- as.character(weightfunction)

  if (binary==F) {

    temp <- data
    temp[temp>0]<-1
    colsum <- colSums(temp)
    colprop <- colsum/nrow(temp)

  }

  else {
    colsum <- colSums(data)
    colprop <- colsum/nrow(data)
  }


  #####
  #####
  ##### Remove columns outside your threshold (and monitor the rows) #####

  if( !(is.null(lower))) {

    data <- data[,which(colsum >= lower)]   # colsum is the sum of nonzero entries

    colsum <- colsum[which(colsum >= lower)]

  }

  if( !(is.null(upper))) {

    data <- data[,which(colsum <= upper)]   # colsum is the sum of nonzero entries

    colsum <- colsum[which(colsum <= upper)]

  }

  rowsum <- rowSums(data)

  if(min(rowsum) <= 0) {  # some rows could lose all nonzero entries when you trim columns

    resp <- readline(prompt="One or more rows has zero weight. \n
                     Make sure that you fix this before continuing. \n
                     Press the ENTER key to continue. \n")

  }




  #####
  #####
  ##### Calculate column-weighted matrix & return #####

  if (sparseinput==F & convertsparse==F) { # if you insist on using a dense matrix

    if (weightfunction == "beta") {   # par1 = alpha, par2 = beta

      x <- seq(0,1, length=1000)
      mode.beta <- max(dbeta(x, shape1=par1, shape2=par2))

      colweights <- dbeta(colprop, shape1=par1, shape2=par2)
      colweights <- colweights/max(mode.beta)    # scale to (0,1) range
      colweights <- sqrt(colweights)
      return(t(t(data)/colweights))

    }

    else if (weightfunction == "step") {    # par1 = min cutoff, par2 = max cutoff

      return(data[,colprop > par1 & colprop < par2])

    }

    else if (weightfunction == "linear") {

      slope1 = 1/mode
      slope2 = -1/(1-mode)

      linweight <- function (density) {
        if (density < mode) { return(slope1*density) }
        else{return(slope2*(density-1) ) }
      }

      colweights <- sapply(colprop, linweight)
      colweights <- sqrt(colweights)
      return(t(t(data)/colweights))

    }

    else if (weightfunction == "IDF") {

      # IDF column weighting = log( N/ 1+density )
      data.idf <- log(nrow(data)/(1 + colsum))
      data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)

      # multiply each column by its IDF weight
      data.tfidf <- crossprod(t(data), data.idf.diag)
      return(data.tfidf)

      # Row normalize
      # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
      # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
      # return(data.tfidf.rn)

    }

    else if (weightfunction == "IDF^2") {

      # IDF column weighting = log( N/ 1+density )
      data.idf <- (log(nrow(data)/(1 + colsum)))^2

      # Multiply each column by its IDF weight
      data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
      data.tfidf <- crossprod(t(data), data.idf.diag)
      return(data.tfidf)

      # Row normalize
      # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
      # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
      # return(data.tfidf.rn)

    }

    else if (weightfunction == "none") {

      return(data)

    }

    else {stop("Pick a valid weight method.")}

  }

  else { # sparse matrix calculations

    if (weightfunction == "beta") {   # par1 = alpha, par2 = beta

      x <- seq(0,1, length=1000)
      mode.beta <- max(dbeta(x, shape1=par1, shape2=par2))

      colweights <- dbeta(colprop, shape1=par1, shape2=par2)
      colweights <- colweights/max(mode.beta)    # scale to (0,1) range
      colweights <- sqrt(colweights)
      return(t(t(data)/colweights))

    }

    else if (weightfunction == "step") {    # par1 = min cutoff, par2 = max cutoff

      return(data[,colprop > par1 & colprop < par2])

    }

    else if (weightfunction == "linear") {

      slope1 = 1/mode
      slope2 = -1/(1-mode)

      linweight <- function (density) {
        if (density < mode) { return(slope1*density) }
        else{return(slope2*(density-1) ) }
      }

      colweights <- sapply(colprop, linweight)
      colweights <- sqrt(colweights)
      return(t(t(data)/colweights))

    }

    else if (weightfunction == "IDF") {

      # IDF column weighting = log( N/ density )
      data.idf <- log(nrow(data)/(colsum))

      # Multiply each column by its IDF weight
      data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
      data.tfidf <- crossprod(t(data), data.idf.diag)
      return(data.tfidf)

      # Row normalize
      # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
      # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
      # return(data.tfidf.rn)

    }

    else if (weightfunction == "IDF^2") {

      # IDF column weighting = log( N/ density )
      data.idf <- (log(nrow(data)/(colsum)))^2

      # Multiply each column by its IDF weight
      data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
      data.tfidf <- crossprod(t(data), data.idf.diag)
      return(data.tfidf)

      # Row normalize
      # data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
      # data.tfidf.rn <- data.tfidf/ rowSums(data.tfidf)
      # return(data.tfidf.rn)

    }

    else if (weightfunction == "none") {

      return(data)

    }

    else {stop("Pick a valid weight method.")}

  }

}
