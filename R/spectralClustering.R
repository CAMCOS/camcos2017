

clustering <- function(Weights, method, k) {
    
    a <- Sys.time()
    
    if (method == "NJW") {
        
        D <- rowSums(Weights)
        
        if( min(D) <= 0) {
            
            resp <- readline(prompt="One of your similarity rows has zero weight. Would you like to set 
                a 1 on the diagonal of the similarity? Type Y or N")
            
            if (resp == "Y" | resp == "y") {
                
                n <- which(D == 0)
                D[n] <- 1 
            }
                
            else { stop("One of your rows has zero weight.") }
            
        }
        
        D <- diag(as.vector(D^-.5))
        Z <- D %*% Weights %*% D
        EZ <- eigen(Z)
        EZ <- EZ$vectors[,1:k]
        U <- t(apply(EZ, 1, function(x) x/sqrt(sum(x^2))))
        NJW <- kmeans(U, centers=k, nstart = 100)
        
        b <- Sys.time()
        time.elapsed <- b - a
        cat("This function took ", time.elapsed, " seconds \n \n")
        
        return(NJW$cluster) }
    
}