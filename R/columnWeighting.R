

# NOTE: this function assumes that multiplying by sqrt(weight) is appropriate

# weightfunction: "beta" "step" "linear"
# par1, par2: beta parameters or indicator function boundaries 
# mode: linear function max (desired density max weight, e.g. 1/k)

colweights <- function (data, weightfunction, par1=NULL, par2=NULL, mode=NULL) {

    temp <- matrix(rep(0, dim(data)[1]*dim(data)[2], nrow=nrow(data)))
    temp[data != 0] <- 1
    colsum <- colSums(temp)
    colprop <- colsum/nrow(temp)
    
    weightfunction <- as.character(weightfunction)
    
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
    
    else {stop("Pick a valid weight method.")}
    
}


