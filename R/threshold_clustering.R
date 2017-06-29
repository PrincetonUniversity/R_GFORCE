
#' Threshold Clustering.
#'
#' Solves the K-means problem using kmeans++ for the initialization and then 
#' runs Lloyd's algorithm.
#' @param X \eqn{n x m} matrix. Each row is treated as a point in \eqn{R^m}.
#' @param K integer. The number of clusters to group the data into.
#' @param R_only logical expression. If \code{R_only == FALSE}, then the included
#' native code implementation will be used. Otherwise, an R implementation is used.
#' @return Returns an object with the components:
#' \describe{
#' \item{\code{clusters}}{a \eqn{n} dimensional integer vector. Entry \eqn{i} to the cluster assignment of the data point given by row \eqn{i} of \code{X}.}
#' \item{\code{centers}}{a \eqn{K x m} numeric matrix. Row \eqn{i} corresponds to the center of cluster \eqn{i}.}
#' \item{\code{num_iters}}{an integer. Number of iterations of Lloyd's Algorithm.}
#' \item{\code{time}}{a numeric. Runtime of Lloyd's Algorithm.}
#' }
#'
#' @examples
#' m <- 10 
#' n <- 10
#' X <- matrix(mvrnorm(m*n,rep(0,m*n),diag(m*n)), nrow = n)
#' km_res <- gforce.kmeans(X,3)
#'
#' @export
gforce.threshold_clustering <- function(X, threshold = 0.7,mode=1){
  res <- NULL
  d <- ncol(X)

  if(mode == 1){
    unclust <- 1:d
    clusters <- list()
    K <- 0
    for(i in 1:d){
      if(i %in% unclust){
        K <- K + 1
        diag_element <- X[i,i]
        i_row <- X[i,]
        same_group <- which(i_row > threshold*diag_element)
        same_group <- intersect(same_group,unclust)
        clusters[[K]] <- same_group
        unclust <- setdiff(unclust,same_group)
      }
    }
  } else {
    unclust <- 1:d
    adj_mat <- matrix(rep(0,d^2),ncol=d)
    for(i in 1:d){
      K <- K + 1
      diag_element <- X[i,i]
      i_row <- X[i,]
      same_group <- which(i_row > threshold*diag_element)
      adj_mat[i,same_group] <- 1
    }
    return(adj_mat)    
  }

  res$K <- K
  res$clusters <- clusters
  return(res)
}
