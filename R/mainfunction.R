#' Main Function to Run Complete Algorithm
#'
#' The accompanying pdf has all the details for this function. Please see this
#' document for details about all function arguments.
#'
#' @param data See description.
#' @param nclust See description.
#' @param sparse See description.
#' @param convertsparse See description.
#' @param save See description.
#' @param return See description.
#' @param preset See description.
#' @param weightfunction See description.
#' @param binary See description.
#' @param lower See description.
#' @param upper See description.
#' @param par1 See description.
#' @param par2 See description.
#' @param mode See description.
#' @param simfunction See description.
#' @param simscale See description.
#' @param rowscaling See description.
#' @param colscaling See description.
#' @param sigma See description.
#' @param centers See description.
#' @param seed See description.
#' @param distance See description.
#' @param clusterfunction See description.
#' @param t See description.
#' @param kmeans.method See description.
#' @param m See description.
#' @param weight.SVD See description.
#' @param SVDdim See description.
#' @param SVDprint See description.
#' @param dim1 See description.
#' @param dim2 See description.
#' @param dim3 See description.
#' @param filepath See description.
#' @param SVDsim See description.
#' @param simdim See description.
#' @param dim1sim See description.
#' @param dim2sim See description.
#' @param dim3sim See description.
#' @param SVDsim.plot See description.
#' @param sim.filepath See description.
#' @param insights See description.
#' @param vocab See description.
#' @param nfeatures See description.
#' @param n.ins See description.
#' @param insight.plot See description.
#' @param insight.filepath See description.
#' @param end.of.arguments See description.
#'
#' @return List.
#' @import clustrd
#' @import corpcor
#' @import dummies
#' @import fclust
#' @import ggplot2
#' @import grid
#' @import irlba
#' @import Matrix
#' @import plyr
#' @import rgl
#' @import RSpectra
#' @import stats
#' @export
mainfunction <- function(data, nclust,

                         # Data input type:
                         sparse=TRUE, convertsparse=TRUE,

                         # Saving and/or Returning the output:
                         save=TRUE, return=TRUE,

                         # Preset argument combinations:
                         preset = 0,

                         # Column weighting argument defaults:
                         weightfunction="IDF", binary=TRUE, lower=2, upper=NULL, par1=NULL, par2=NULL, mode=NULL,

                         # Similarity function argument defaults:
                         simfunction = "cosine", simscale=NULL,
                         rowscaling = NULL, colscaling = NULL, sigma = NULL, centers = NULL, seed = NULL, distance = NULL,

                         # Clustering function argument defaults:
                         clusterfunction="DiffusionMap", t=.5, kmeans.method="kmeans", m=NULL,

                         # SVD options for weighted (IDF) data:
                         weight.SVD=FALSE, SVDdim=200, SVDprint=FALSE,
                         dim1=1, dim2=2, dim3=3, filepath=NULL, # for plotting 3D graphs
                         # specify filepath if you want to save as PDF

                         # SVD options for similarity (cosine) matrix:
                         SVDsim=TRUE, simdim=3, dim1sim=1, dim2sim=2, dim3sim=3, SVDsim.plot=FALSE, sim.filepath=NULL,

                         # Cluster Insights:
                         insights=FALSE, vocab=NULL, nfeatures=20, n.ins=NULL, insight.plot = TRUE, insight.filepath = NULL,
                         # insight.filepath should be a folder, not a filename, if you want to store multiple files


                         end.of.arguments=NULL) # end of arguments (for neatness)

{

  # newfolder <- gsub(":", "_", paste ("~/",Sys.time(),sep="_"))
  #
  # dir.create(newfolder)
  #
  # setwd(newfolder)

  require(Matrix)


  #####
  #####
  ##### PRESET ARGUMENTS #####

  if (preset==1) {

    weight.SVD = TRUE
    weightfunction = "IDF^2"
  }

  if (preset==2) {

    kmeans.method="poly.fuzzy"
  }

  if (preset==3) {

    clusterfunction="NJW"
    kmeans.method="poly.fuzzy"
  }

  if (preset==4) {

    weight.SVD = TRUE
    weightfunction = "IDF^2"
    kmeans.method="poly.fuzzy"
  }

  if (preset==5) {

    weight.SVD = TRUE
    weightfunction = "IDF^2"
    clusterfunction="NJW"
    kmeans.method="poly.fuzzy"
  }


  #####
  #####
  ##### DEFINE OUR FUNCTIONS #####

  colweights <- function (data, weightfunction, sparseinput,
                          par1=NULL, par2=NULL, mode=NULL,
                          binary=TRUE, convertsparse=TRUE,
                          lower=2, upper=NULL) {

    #####
    ##### construct Matrix object for use with "Matrix" package #####

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

      # resp <- readline(prompt="One or more rows has zero weight. \n
      #     Make sure that you fix this before continuing. \n
      #     Press the ENTER key to continue. \n")

      badrows <- which(rowsum<=0)

      data <- data[-badrows,]

      cat(length(badrows), " rows have zero weight, and will be removed.")



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


  similarity <- function(data, method, rowscaling = NULL, colscaling = NULL,
                         sigma = NULL, centers = NULL, seed = NULL, distance = NULL,
                         sparse = T, simscale) {

    #####
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
    ##### Distance -> Gaussian similarity (if applicable) #####

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
    ##### Compute the NxN similarity matrix and return #####

    ###
    ### dense matrix:
    ###

    if (sparse == F) {

      if (is.null(centers)) {

        if (method == "correlation") {

          centeredcolumns <- t(t(data)-colMeans(data)) # center the data by column

          rowvar <- rowSums(centeredcolumns^2) # store the variances for each row (where colMeans=0)

          Cov.matrix <- tcrossprod(centeredcolumns) # calculate the covariance matrix (dot product all rows)

          # corr = cov(x,y) / sqrt(var(x)var(y)
          Corr.true <- Cov.matrix/sqrt(rowvar)
          Similarity <- t(t(Corr.true)/sqrt(rowvar))

          # # scale the matrix to (0,1), excluding the diagonal
          # diag(Similarity) <- rep(0,nrow(Similarity))
          # Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
          # diag(Similarity) <- rep(0, nrow(Similarity))

          # Set any negative similarities equal to zero

          Similarity[Similarity<0] <- 0

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

          # # scale the matrix to (0,1), excluding the diagonal
          # diag(Similarity) <- rep(0,nrow(Similarity))
          # Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
          # diag(Similarity) <- rep(0, nrow(Similarity))

          # Set any negative similarities equal to zero

          Similarity[Similarity<0] <- 0

        }

        else if (method == "cosine") {

          rowlength <- rowSums(data^2)
          Dot.prods <- tcrossprod(data)
          Cosines <- Dot.prods/sqrt(rowlength)
          Similarity <- t(t(Cosines)/sqrt(rowlength))
          # diag(Similarity) <- rep(0,nrow(Similarity))
          # Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
          # diag(Similarity) <- rep(0, nrow(Similarity))

          # Set any negative similarities equal to zero

          Similarity[Similarity<0] <- 0

        }

        else if (method == "dotproduct") {

          Similarity <- tcrossprod(data)
          # diag(Similarity) <- rep(0,nrow(Similarity))
          # Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
          # diag(Similarity) <- rep(0, nrow(Similarity))

          # Set any negative similarities equal to zero

          Similarity[Similarity<0] <- 0

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

    ###
    ### sparse matrix:
    ###

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

          if(is.null(simscale)) {

            if (weight.SVD==TRUE) {simscale="negative"}

            if (weight.SVD==FALSE) {simscale="0-1"}
          }

          if (simscale == "negative") {

            # Set any negative similarities equal to zero

            Similarity[Similarity<0] <- 0
            diag(Similarity) <- rep(0, nrow(Similarity))
          }

          else if (simscale == "0-1") {

            # rescale to (0,1) interval

            diag(Similarity) <- rep(0,nrow(Similarity))
            Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
            diag(Similarity) <- rep(0, nrow(Similarity))
          }

        }

        else if (method == "dotproduct") {

          Similarity <- tcrossprod(data)

          if(is.null(simscale)) {

            if (weight.SVD==TRUE) {simscale="negative"}

            if (weight.SVD==FALSE) {simscale="0-1"}
          }

          if (simscale == "negative") { # Set negative similarities to zero

            Similarity[Similarity<0] <- 0
          }

          else if (simscale == "0-1") { # rescale to (0,1) interval

            diag(Similarity) <- rep(0,nrow(Similarity))
            Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
            diag(Similarity) <- rep(0, nrow(Similarity))
          }

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

    return(Similarity)

  }


  clustering <- function(Weights, k, method, t=NULL, sparse=T,
                         kmeans.method="kmeans", m=NULL) {

    if (method == "NJW") {

      D <- rowSums(Weights)   # Degrees matrix #####

      ###
      ### Check: if row has zero similarity, problems arise
      ###

      if( min(D) <= 0) {

        resp <- readline(prompt="One of your similarity rows has zero weight. Would you like to set
                         a 1 on the diagonal of the similarity? Type Y or N \n")
        if (resp == "Y" | resp == "y") {
          n <- which(D == 0)
          D[n] <- 1
        }

        else { stop("One of your rows has zero weight.") }
      }

      #####
      ##### Subspace Projection #####
      #####

      D <- Diagonal(n=nrow(Weights),(D^-.5))
      Z <- D %*% Weights %*% D


      if (kmeans.method=="RKM") {

        ###
        ### Define the RKM function
        ###

        RKM <- function (data, nclus, ndim, alpha = NULL, method = "RKM", center = TRUE,
                         scale = TRUE, rotation = "none", nstart = 10, smartStart = NULL,
                         seed = 1234) {

          require(ggplot2)
          require(dummies)
          require(grid)
          require(corpcor)

          ssq = function(a) {
            t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
          }

          if (is.null(alpha) == TRUE) {
            if (method == "RKM") {
              alpha = 0.5
            }
            else if (method == "FKM") {
              alpha = 0
            }
          }
          odata = data
          data = scale(data, center = center, scale = scale)
          # data = data.matrix(data)
          n = dim(data)[1]
          m = dim(data)[2]
          conv = 1e-06
          func = {
          }
          index = {
          }
          AA = {
          }
          FF = {
          }
          YY = {
          }
          UU = {
          }

          require(irlba)

          for (run in c(1:nstart)) {
            if (is.null(smartStart)) {
              myseed = seed + run
              set.seed(myseed)
              randVec = matrix(ceiling(runif(n) * nclus), n, 1)
            }
            else {
              randVec = smartStart
            }
            U = dummy(randVec)
            P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)


            # A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) *
            #         diag(n)) %*% data)$vectors
            # A = A[, 1:ndim]


            testobj <- t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * diag(n)) %*% data

            A <- partial_eigen(x=testobj, n = ndim, symmetric = TRUE)$vectors


            G = data %*% A
            Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
            f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * ssq(data %*%
                                                                     A - U %*% Y)
            f = as.numeric(f)
            fold = f + 2 * conv * f
            iter = 0
            while (f < fold - conv * f) {
              fold = f
              iter = iter + 1
              outK = try(kmeans(G, centers = Y, nstart = 100),
                         silent = T)
              if (is.list(outK) == FALSE) {
                outK = EmptyKmeans(G, centers = Y)
              }
              v = as.factor(outK$cluster)
              U = diag(nlevels(v))[v, ]
              P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)

              # A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 *
              #         alpha) * diag(n)) %*% data)$vectors
              # A = A[, c(1:ndim)]

              testobj <- t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * diag(n)) %*% data

              A <- partial_eigen(x=testobj, n = ndim, symmetric = TRUE)$vectors


              G = data %*% A
              Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
              f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) *
                ssq(data %*% A - U %*% Y)
            }
            func[run] = f
            FF[[run]] = G
            AA[[run]] = A
            YY[[run]] = Y
            UU[[run]] = U
            cat("Just finished iteration ", run, "\n")
          }
          mi = which.min(func)
          U = UU[[mi]]
          cluID = apply(U, 1, which.max)
          csize = round((table(cluID)/sum(table(cluID))) * 100, digits = 2)
          aa = sort(csize, decreasing = TRUE)
          require(plyr)
          cluID = mapvalues(cluID, from = as.integer(names(aa)), to = as.integer(names(table(cluID))))
          centroid = YY[[mi]]
          centroid = centroid[as.integer(names(aa)), ]
          if (rotation == "varimax") {
            require(stats)
            AA[[mi]] = varimax(AA[[mi]])$loadings
            FF[[mi]] = data %*% AA[[mi]]
            centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% FF[[mi]]
            centroid = centroid[as.integer(names(aa)), ]
          }
          else if (rotation == "promax") {
            AA[[mi]] = promax(AA[[mi]])$loadings[1:m, 1:ndim]
            FF[[mi]] = data %*% AA[[mi]]
            centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% FF[[mi]]
            centroid = centroid[as.integer(names(aa)), ]
          }
          out = list()
          mi = which.min(func)
          out$obscoord = FF[[mi]]
          rownames(out$obscoord) = rownames(data)
          out$attcoord = data.matrix(AA[[mi]])
          rownames(out$attcoord) = colnames(data)
          out$centroid = centroid
          names(cluID) = rownames(data)
          out$cluID = cluID
          out$criterion = func[mi]
          out$csize = round((table(cluID)/sum(table(cluID))) * 100,
                            digits = 1)
          out$odata = odata
          out$scale = scale
          out$center = center
          out$nstart = nstart
          class(out) = "cluspca"
          return(out)
        }


        ###
        ### Run RKM on the normalized Z matrix
        ###

        cluster.out = RKM(Z, nclus=k, ndim=k, method = "RKM", rotation = "varimax", nstart=10)


        return(cluster.out)
      }


      # RSpectra is efficient for dense matrices


      if (k=="eigenvalue") {

        require(irlba)

        k <- 20

        EZ <- partial_eigen(x=Z, n = k, symmetric = TRUE)$value

        print(EZ)

        k <- readline(prompt="\n Here is a list of the first 20 eigenvalues.
                      Pick whichever eigenvalue appears best. \n")

        k <- as.numeric(k)

      }


      if (sparse==F) {
        require(RSpectra)
        EZ <- eigs_sym(Z, k+1, 'LM')$vector
        EZ <- EZ[,1:k]
      }


      # irlba is efficient and accurate for sparse matrices

      else {
        require(irlba)
        EZ <- partial_eigen(x=Z, n = k+1, symmetric = TRUE)$vectors
        EZ <- EZ[,1:k]
      }


      # U is the L2-normalized eigenspace

      U <- EZ/sqrt(rowSums(EZ^2))


      #####
      ##### k-Mmeans in this normalized eigenspace:
      #####
      {
        # Regular k-means
        if (kmeans.method=="kmeans") {
          cluster.out <- kmeans(U, centers=k, nstart = 100)
        }

        # Fuzzy k-means
        else if (kmeans.method=="fuzzy") {
          require(fclust)
          if (is.null(m)) {
            m <- 2
          }
          cluster.out <- FKM(X=U,k=k, m=m, RS=10)
        }

        # Polynomial fuzzy k-means
        else if (kmeans.method=="poly.fuzzy") {
          require(fclust)
          if (is.null(m)) {
            m <- .5
          }
          cluster.out <- FKM.pf(X=U,k=k, b=m, RS=10)
        }

        else {stop("Pick a valid k-means method.")}
      }

  }


    else if (method == "Ncut") {

      n <- nrow(Weights)
      dvec_inv = 1/sqrt(rowSums(Weights))
      #W_tilde = Matrix(rep(dvec_inv,n), ncol=n) * Weights * t(Matrix(rep(dvec_inv,n),ncol=n))
      W_tilde = Diagonal(n,dvec_inv) %*% Weights %*% Diagonal(n,dvec_inv)
      W_tilde = (W_tilde+t(W_tilde))/2

      # diag(dvec_inv) %*% Weights %*% diag(dvec_inv) ?
      # why the average part?

      if (sparse==F) {
        require(RSpectra)
        EZ <- eigs_sym(W_tilde, k, 'LM')$vector }
      else {
        require(irlba)
        EZ <- partial_eigen(x=W_tilde, n = k, symmetric = TRUE)$vectors
      }

      V <- EZ
      V = matrix(rep(dvec_inv,k-1), ncol = k-1) * V[,2:k]
      V = V / (matrix(rep(sqrt(rowSums(V^2)),k-1),ncol=k-1))


      if (kmeans.method=="kmeans") {
        cluster.out <- kmeans(V, centers=k, nstart = 100)
      }

      else if (kmeans.method=="fuzzy") {
        require(fclust)
        if (is.null(m)) {
          m <- 2
        }
        cluster.out <- FKM(X=V,k=k, m=m, RS=10)
      }

      else if (kmeans.method=="poly.fuzzy") {
        require(fclust)
        if (is.null(m)) {
          m <- .5
        }
        cluster.out <- FKM.pf(X=V,k=k, b=m, RS=10)
      }

      else {stop("Pick a valid k-means method.")}

    }


    else if (method == 'DiffusionMap'){

      if (is.null(t)) {

        stop("Specify a t value.") }

      require(RSpectra)

      n <- nrow(Weights)
      dvec_inv = 1/sqrt(rowSums(Weights))
      #W_tilde = matrix(rep(dvec_inv,n), ncol=n) * Weights * t(matrix(rep(dvec_inv,n),ncol=n))
      W_tilde = Diagonal(n,dvec_inv) %*% Weights %*% Diagonal(n,dvec_inv)
      W_tilde = (W_tilde+t(W_tilde))/2

      if (kmeans.method=="RKM") {

        ###
        ### Define the RKM function
        ###

        RKM <- function (data, nclus, ndim, alpha = NULL, method = "RKM", center = TRUE,
                         scale = TRUE, rotation = "none", nstart = 10, smartStart = NULL,
                         seed = 1234) {

          require(ggplot2)
          require(dummies)
          require(grid)
          require(corpcor)

          ssq = function(a) {
            t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
          }

          if (is.null(alpha) == TRUE) {
            if (method == "RKM") {
              alpha = 0.5
            }
            else if (method == "FKM") {
              alpha = 0
            }
          }
          odata = data
          data = scale(data, center = center, scale = scale)
          # data = data.matrix(data)
          n = dim(data)[1]
          m = dim(data)[2]
          conv = 1e-06
          func = {
          }
          index = {
          }
          AA = {
          }
          FF = {
          }
          YY = {
          }
          UU = {
          }

          require(irlba)

          for (run in c(1:nstart)) {
            if (is.null(smartStart)) {
              myseed = seed + run
              set.seed(myseed)
              randVec = matrix(ceiling(runif(n) * nclus), n, 1)
            }
            else {
              randVec = smartStart
            }
            U = dummy(randVec)
            P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)


            # A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) *
            #         diag(n)) %*% data)$vectors
            # A = A[, 1:ndim]


            testobj <- t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * diag(n)) %*% data

            A <- partial_eigen(x=testobj, n = ndim, symmetric = TRUE)$vectors


            G = data %*% A
            Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
            f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) * ssq(data %*%
                                                                     A - U %*% Y)
            f = as.numeric(f)
            fold = f + 2 * conv * f
            iter = 0
            while (f < fold - conv * f) {
              fold = f
              iter = iter + 1
              outK = try(kmeans(G, centers = Y, nstart = 100),
                         silent = T)
              if (is.list(outK) == FALSE) {
                outK = EmptyKmeans(G, centers = Y)
              }
              v = as.factor(outK$cluster)
              U = diag(nlevels(v))[v, ]
              P = U %*% pseudoinverse(t(U) %*% U) %*% t(U)

              # A = eigen(t(data) %*% ((1 - alpha) * P - (1 - 2 *
              #         alpha) * diag(n)) %*% data)$vectors
              # A = A[, c(1:ndim)]

              testobj <- t(data) %*% ((1 - alpha) * P - (1 - 2 * alpha) * diag(n)) %*% data

              A <- partial_eigen(x=testobj, n = ndim, symmetric = TRUE)$vectors


              G = data %*% A
              Y = pseudoinverse(t(U) %*% U) %*% t(U) %*% G
              f = alpha * ssq(data - G %*% t(A)) + (1 - alpha) *
                ssq(data %*% A - U %*% Y)
            }
            func[run] = f
            FF[[run]] = G
            AA[[run]] = A
            YY[[run]] = Y
            UU[[run]] = U
            cat("Just finished iteration ", run, "\n")
          }
          mi = which.min(func)
          U = UU[[mi]]
          cluID = apply(U, 1, which.max)
          csize = round((table(cluID)/sum(table(cluID))) * 100, digits = 2)
          aa = sort(csize, decreasing = TRUE)
          require(plyr)
          cluID = mapvalues(cluID, from = as.integer(names(aa)), to = as.integer(names(table(cluID))))
          centroid = YY[[mi]]
          centroid = centroid[as.integer(names(aa)), ]
          if (rotation == "varimax") {
            require(stats)
            AA[[mi]] = varimax(AA[[mi]])$loadings
            FF[[mi]] = data %*% AA[[mi]]
            centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% FF[[mi]]
            centroid = centroid[as.integer(names(aa)), ]
          }
          else if (rotation == "promax") {
            AA[[mi]] = promax(AA[[mi]])$loadings[1:m, 1:ndim]
            FF[[mi]] = data %*% AA[[mi]]
            centroid = pseudoinverse(t(U) %*% U) %*% t(U) %*% FF[[mi]]
            centroid = centroid[as.integer(names(aa)), ]
          }
          out = list()
          mi = which.min(func)
          out$obscoord = FF[[mi]]
          rownames(out$obscoord) = rownames(data)
          out$attcoord = data.matrix(AA[[mi]])
          rownames(out$attcoord) = colnames(data)
          out$centroid = centroid
          names(cluID) = rownames(data)
          out$cluID = cluID
          out$criterion = func[mi]
          out$csize = round((table(cluID)/sum(table(cluID))) * 100,
                            digits = 1)
          out$odata = odata
          out$scale = scale
          out$center = center
          out$nstart = nstart
          class(out) = "cluspca"
          return(out)
        }


        ###
        ### Run RKM on the normalized W_tilde matrix
        ###

        cluster.out = RKM(W_tilde, nclus=k, ndim=k+1, method = "RKM", rotation = "varimax", nstart=10)

        return(cluster.out)
      }

      require(irlba)

      if (k=="eigenvalue") {

        k <- 20

        EZ <- partial_eigen(x=W_tilde, n = k, symmetric = TRUE)$value

        print(EZ)

        k <- readline(prompt="\n Here is a list of the first 20 eigenvalues.
                      Pick whichever eigenvalue appears best. \n")

        k <- as.numeric(k)

      }

      EV <- eigs_sym(W_tilde, k+1, 'LM')
      V <- EV$vector
      lambda <- EV$value

      V_inv = 1/sqrt((rowSums(V[,2:(k+1)]^2)))
      V <- matrix(rep(V_inv,k), ncol=k) * V[,2:(k+1)]
      V = matrix(rep(dvec_inv,k), ncol = k)  * V

      V <-  abs(matrix(rep(lambda[2:(k+1)], each=n), ncol=k))^(t)* V

      # V = (matrix(rep(lambda[2:(k)], each=n), ncol=k-1)^t )* V

      # if( (t%%1) != 0){
      #     V <-  abs(matrix(rep(lambda[1:(k)], each=n), ncol=k)^(t))* V
      # }
      # else {
      #     V <-  abs(matrix(rep(lambda[1:(k)], each=n), ncol=k)^(t))* V
      # }

      # run kmeans in eigenspace:

      if (kmeans.method=="kmeans") {
        cluster.out <- kmeans(V, centers=k, nstart = 100)
      }

      else if (kmeans.method=="fuzzy") {
        require(fclust)
        if (is.null(m)) {
          m <- 2
        }
        cluster.out <- FKM(X=V,k=k, m=m, RS=10)
      }

      else if (kmeans.method=="poly.fuzzy") {
        require(fclust)
        if (is.null(m)) {
          m <- .5
        }
        cluster.out <- FKM.pf(X=V,k=k, b=m, RS=10)
      }

      else if (kmeans.method=="RKM") {
        require(clustrd)
        cluster.out = cluspca(V, nclus=k, ndim=k, method = "RKM", rotation = "varimax", nstart=10)
      }

    }


    else {stop("Pick a valid clustering method.") }

  }


  #####
  #####
  ##### RUN OUR FUNCTIONS #####


  ### Store a copy of the data for summary statistics later ###

  copydata <- colweights(data, weightfunction="none", sparseinput=sparse, binary=T)


  ### Column weighting ###

  col.args = list(data=data, weightfunction=weightfunction, sparseinput=sparse,
                  par1=par1, par2=par2, mode=mode,
                  binary=binary, convertsparse=convertsparse,
                  lower=lower, upper=upper)


  weighteddata <- do.call(colweights, col.args)


  cat("Column weighting is finished. \n")


  ### if we want to get insights later, we need to store a copy of the data

  if (insights==TRUE & weight.SVD==TRUE){weighteddata2 <- weighteddata}


  ### change "sparse" to true if you converted to sparse in colweights step:

  if (convertsparse==T) {sparse=T}


  ### Do you want to calculate SVD on the weighted data? ###

  svd.data = NULL

  if (weight.SVD==TRUE) {

    require(irlba)
    all.svd <- irlba(weighteddata, SVDdim)
    svd.data <- weighteddata %*% all.svd$v

    if (SVDprint==TRUE) {   # print the SVD results

      require(rgl)
      svdoutput <- plot3d(svd.data[,dim1], svd.data[,dim2], svd.data[,dim3], col="blue")

      if (!(is.null(filepath))) {

        snapshot3d(filepath)
      }
    }

    weighteddata <- svd.data

    sparse = F

    convertsparse = F

    cat("SVD is finished. \n")

  }


  ### Similarity Matrix ###

  sim.args = list(data=weighteddata, method=simfunction,
                  rowscaling = rowscaling, colscaling = colscaling,
                  sigma = sigma, centers = centers, seed = seed,
                  distance = distance, sparse = sparse, simscale=simscale)

  simdata <- do.call(similarity, sim.args)

  diag(simdata) <- 0


  cat("Similarity matrix is finished. \n")


  ### Do you want to plot the SVD of the Similarity matrix? ###

  simdata.svd = NULL

  if (SVDsim==TRUE) {

    require(irlba)
    sim.svd <- irlba(simdata, simdim)
    simdata.svd <- simdata %*% sim.svd$v

    if (SVDsim.plot==T) {
      require(rgl)
      svdoutput <- plot3d(simdata.svd[,dim1sim], simdata.svd[,dim2sim], simdata.svd[,dim3sim], col="blue")
    }

    if (!(is.null(sim.filepath))) {

      snapshot3d(sim.filepath)
    }

    cat("Similarity SVD is finished. \n")


  }


  ### Clustering step ###

  clust.args = list(Weights=simdata, k=nclust, method=clusterfunction,
                    t=t, sparse=sparse, kmeans.method=kmeans.method, m=m)


  clusterdata <- do.call(clustering, clust.args)


  cat("Clustering is finished. \n")


  #####
  #####
  ##### Cluster Insights #####

  if (insights==TRUE & weight.SVD==TRUE){weighteddata <- weighteddata2}

  insight.output = NULL

  if (insights==TRUE) {

    getInsights <- function(cluster, vocab, n,
                            plot, file){

      require(RSpectra)
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

      return(list(magnitude = v1_topn, keywords = as.character(v1_top.words)))

    }


    if (is.null(n.ins)) { n.ins <- nclust }

    insight.output = list()

    for (clustID in 1:n.ins) {

      insight.output[[n.ins]] <- getInsights(cluster=weighteddata[clusterdata$clus==clustID,], vocab=vocab,
                                             n=nfeatures, plot=T, file=insight.filepath)

    }

    cat("Cluster insights is finished. \n")

  }




  #####
  #####
  ##### Gather the Output #####

  # Gather any statistics of interest:

  summarystats <- list(
    dimensions=dim(copydata),
    colsum=colSums(copydata),
    rowsum=rowSums(copydata),
    density=nnzero(copydata),
    max(copydata)
  )


  # Remove any data from the output:

  col.args <- col.args[-1]
  sim.args <- sim.args[-1]
  clust.args <- clust.args[-1]


  # Output a list of relevant objects:

  output <- list(col.args=col.args, sim.args=sim.args, clust.args=clust.args,
                 summarystats=summarystats, clusterdata=clusterdata,
                 svd.data=svd.data, simdata.svd=simdata.svd, insight.output=insight.output
  )#[which(c(T,T,T,T,T, weight.SVD, SVDsim, insights))]

  if (save==T) {save(output, file = gsub(":", "_", paste ("~/",Sys.time(),sep="_")) )}

  if (return==T) {return(output)}

  # test <- list(w,x,y,z)[which(c(a,b,c,d))]


}
