#' Checking cluster accuracy
#'
#' in the results table, true clusters are the rows, guessed clusters are the columns
#'
#' @param cluster
#' @param Labels
#' @param k
#'
#' @return
#' @export
#'
#' @examples
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


# clean up the assignment matrix

assignment <- cbind( seq(1,k), assignment[,1])
colnames(assignment) <- c("True Cluster", "Assigned Cluster")


return(list(results, assignment, percent.score))

}

# example of a results matrix that FAILS clustercheck:
#
# results <- matrix( c(15,31,4,452,5,0,1,384,1,20,1,3,34,141,56,53,289,357,0,2,282,6,160,5,527,32,48,51,9,6,4,2,203,0,0,5), nrow=6)
#
# cluster assignment of "true cluster #3" is split too evenly



