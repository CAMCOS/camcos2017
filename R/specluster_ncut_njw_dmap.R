# Input: X : distance matrix;
#		 k : number of clusters;
#		 sigmaValue: sigma for Gaussian simmilarity
#		 t: diffusionmap steps (default is 0)
#		 clustermethod: "NJW", "Ncut", "DiffusionMap"
		 
# Output: clustering results: output$clusters


specluster_ncut_njw_dmap <- function(X, k, sigmaValue, t=0, clusterMethod) {
    X <- as.matrix(X)
    n = nrow(X)
    sigma2 = 2*sigmaValue^2
    
    W = exp(-X^2 / sigma2)
    diag(W) = 0
    #image(W, main= 'weight matrix')
    
    dvec_inv = 1/sqrt(rowSums(W))
    W_tilde = matrix(rep(dvec_inv,n), ncol=n)*W*t(matrix(rep(dvec_inv,n),ncol=n))
    W_tilde = (W_tilde+t(W_tilde))/2
    
    V <- eigs_sym(W_tilde, k, 'LM')$vector
    lambda <- eigs_sym(W_tilde, k, 'LM')$value
    
    if (clusterMethod == 'Ncut'){
        V = matrix(rep(dvec_inv,k-1), ncol = k-1)  * V[,2:k]
        #MATLAB: V = repmat(dvec_inv, 1, k-1).%*% V[,2:k]
    }
    
    if (clusterMethod == 'DiffusionMap'){
        V_inv = 1/sqrt((rowSums(V[,2:k]^2)))
        V <- matrix(rep(V_inv,k-1), ncol=k-1) * V[,2:k]
        V = matrix(rep(dvec_inv,k-1), ncol = k-1)  * V
        V = (matrix(rep(lambda[2:k], each=n), ncol=k-1)^t )* V
        #MATLAB: V = repmat(dvec_inv, 1, k-1).%*% V[,2:k]
    }
    
    if (clusterMethod == 'NJW'){
        V = V / (matrix(rep(sqrt(rowSums(V^2)),k),ncol=k))
        # MATLAB: V = V./repmat(sqrt(sum(V.^2,2)), 1,k)
    }
    labels = kmeans(V, k, 10)
    labels
}

