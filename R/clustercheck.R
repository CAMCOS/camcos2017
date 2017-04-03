#'Check accuracy of cluster assignments
#'
#'Calculate accuracy of the clustering algorithm by comparing algorithm
#'generated cluster labels to the ground truth cluster labels.
#'
#'In the results matrix, ground truh cluster labels are rows and the algorithm
#'generated cluster labels are columns.
#'@param cluster Integer vector of algorithm generated cluster labels.
#'@param Labels Integer vector of ground truth cluster labels.
#'@param k Integer indicating the number of clusters.
#'
#'@return List with two elements:
#'  \enumerate{
#'    \item Matrix of counts indicating how the ground truth labels are
#'    distributed across the algorithm generated cluster labels.
#'    \item Numeric. Percentage correctly assigned.
#'  }
#'
#' @examples
#' data(cluster6)
#' # Perfect cluster assignment
#' clustercheck(cluster6_labels, cluster6_labels, 6)
#'
#' # Random cluster assignment
#' n <- length(cluster6_labels)
#' random_labels <- sample(1:6, size = n, replace = TRUE)
#' clustercheck(random_labels, cluster6_labels, 6)
clustercheck <- function (cluster, Labels, k) {


  comp <- cbind(cluster, Labels) # nx2 matrix of cluster assignments and true group membership

  # results matrix

  results <- matrix( rep(0,k^2), nrow = k)


  # results matrix: rows = true group, columns = assigned cluster

  for (row in 1:nrow(comp)){

      results[comp[row,2], comp[row,1]] <- results[comp[row,2], comp[row,1]] + 1

  }


  # calculate percentage distribution of each true group

  percent.cluster <- t(apply(results, 1, function(x) x/sum(x)))

  assignment <- matrix( rep(0,k), nrow = k)


  # assign the most-correct true group to that cluster, and continue

  for (index in 1:k) {

      test <- which(percent.cluster == max(percent.cluster), arr.ind = TRUE)

      assignment[ test[1],1 ] <- test[2]

      percent.cluster[test[1],] <- 0   # remove that row (true group) from consideration

      percent.cluster[,test[2]] <- 0   # remove that column (cluster) from consideration

  }


  # calculate percentage assigned

  correct <- 0
  for (i in 1:nrow(assignment)){

      correct <- correct + results[i,assignment[i,1]]

  }

  percent.score <- correct / sum(results)

  return(list(results,percent.score))

}

# example of a results matrix that FAILS clustercheck:
#
# results <- matrix( c(15,31,4,452,5,0,1,384,1,20,1,3,34,141,56,53,289,357,0,2,282,6,160,5,527,32,48,51,9,6,4,2,203,0,0,5), nrow=6)
#
# cluster assignment of "true cluster #3" is split too evenly



