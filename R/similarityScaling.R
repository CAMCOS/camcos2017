

### Mandatory:
# data matrix
# method: "correlation" "cosine" "dotproduct" "Gaussian"

### Optional:
# colscaling: "standardize" (center and divide by variance)
# rowscaling: "L2" "L1"
# centers: % of columns to be sampled, if desired (e.g. .25)
# seed: for repeatability
# distance: for use with Gaussian Kernel; options = JSdivergence, L1, L2


similarity <- function(data, method, rowscaling = NULL, colscaling = NULL, 
    sigma = NULL, centers = NULL, seed = NULL, distance = NULL) {
    
    a <- Sys.time()
    
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
        
        else { Similarity <- exp(-1 * dist / (2*sigma)) }
            
    }      
    
    #####
    #####
    ##### Full data similarity matrix #####
    
    if (is.null(centers)) {
    
        if (method == "correlation") {
            
            library(proxy)
            Similarity <- simil(data, method = "correlation")
            Similarity <- as.matrix(Similarity)
            diag(Similarity) <- rep(1,nrow(Similarity)) 
            Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1]) }
        
        else if (method == "cosine") {
            
            library(proxy)
            Similarity <- simil(data, method = "cosine")
            Similarity <- as.matrix(Similarity)
            diag(Similarity) <- rep(1,nrow(Similarity)) }
        
        else if (method == "dotproduct") {
            
            Similarity <- data %*% t(data) }
        
    }
    
    
    #####
    #####
    ##### Random centers similarity matrix #####
    
    else {
        
        if (is.numeric(seed)) { set.seed(seed) }
    
        kcenters <- centers*nrow(data)
        rcenters <- as.matrix(data[sample(nrow(data),kcenters,replace = FALSE), ])
        
        # rcenters is rxN matrix, r = kcenters % of rows, sampled randomly
        
        if (method == "correlation") {
            
            library(proxy)
            Similarity <- simil(data, rcenters, method = "correlation")
            Similarity <- as.matrix(Similarity)
            Similarity <- (Similarity - min(Similarity))/(range(Similarity)[2] - range(Similarity)[1])
            Similarity <- Similarity %*% t(Similarity)  }
        
        else if (method == "cosine") {
            
            library(proxy)
            Similarity <- simil(data, rcenters, method = "cosine")
            Similarity <- as.matrix(Similarity)
            Similarity <- Similarity %*% t(Similarity)  }
        
        else if (method == "dotproduct") {
            
            Similarity <- data %*% t(rcenters) 
            Similarity <- Similarity %*% t(Similarity)  }
        
    }
    
    #####
    #####
    
    b <- Sys.time()
    time.elapsed <- b - a
    
    print(time.elapsed)
    
    return(Similarity)
    
}
    

