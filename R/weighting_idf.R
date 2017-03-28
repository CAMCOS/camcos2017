## Input: data:
## Ouput: row-normalized, column idf-weighted data

weighting_idf <- function(data) {
    data.m <- as.matrix(data)
    ### IDF: Inverse document frequency #####
    data.idf <- log(nrow(data.m)/(1 + colSums(data.m != 0)))
    data.idf.diag <- Diagonal(n = length(data.idf), x=data.idf)
    ## weighted matrix
    data.tfidf <- crossprod(t(data.m), data.idf.diag)
    ## row normalize
    data.tfidf.rn <- data.tfidf/ sqrt(rowSums(data.tfidf^2))
}
