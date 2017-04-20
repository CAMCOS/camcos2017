

# NOTE: this function assumes that multiplying by sqrt(weight) is appropriate

# weightfunction: "beta" "step" "linear" "IDF"
# par1, par2: beta parameters or step function boundaries
# mode: linear function max (desired density max weight, e.g. 1/k)

colweights <- function (data, weightfunction, par1=NULL, par2=NULL,
    mode=NULL, binary=TRUE, sparse=TRUE) {

    require(Matrix)

    # construct Matrix object for use with Matrix package #####
    ### converts matrix and data frame object to Matrix

    if (sparse==T) {    # given a sparse matrix

        if (is.matrix(data)) {
            if (binary==T) {
                data <- sparseMatrix(i=data[,1], j=data[,2])
            }

            else {
                data <- sparseMatrix(i=data[,1], j=data[,2], x=data[,3])
            }
        }

        else if (is.data.frame(data)) {

            data <- as.matrix(data)
            data <- Matrix(data)

            if (binary==T) {
                data <- sparseMatrix(i=data[,1], j=data[,2])
            }

            else {
                data <- sparseMatrix(i=data[,1], j=data[,2], x=data[,3])
            }
        }
    }

    else{   # given a dense matrix

        if (is.matrix(data)) {
            data <- Matrix(data)
        }

        else if (is.data.frame(data)) {
            data <- as.matrix(data)
            data <- Matrix(data)
        }

    }



    # find density proportion of each column #####

    weightfunction <- as.character(weightfunction)

    if( !(weightfunction=="IDF")) { # idf doesn't use column proportions

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



    }


    # calculate similarity matrix #####

    if (sparse==F) {

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

            ### IDF: Inverse document frequency #####
            data.idf <- log(nrow(data)/(1 + colSums(data)))
            data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
            ## weighted matrix
            data.tfidf <- crossprod(t(data), data.idf.diag)
            ## row normalize
            data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
            return(data.tfidf.rn)

        }

        else {stop("Pick a valid weight method.")}

    }

    else {

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

            ### IDF: Inverse document frequency #####
            data.idf <- log(nrow(data)/(1 + colSums(data)))
            data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
            ## weighted matrix
            data.tfidf <- crossprod(t(data), data.idf.diag)
            ## row normalize
            data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
            return(data.tfidf.rn)

        }

        else {stop("Pick a valid weight method.")}

    }

}


