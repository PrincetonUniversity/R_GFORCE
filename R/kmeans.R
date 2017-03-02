# Interface to C implementation of K-means++ and Lloyd's Algorithm

# centers are the columns

#' K-means Clustering.
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
#' }
#'
#' @examples
#' m <- 10 
#' n <- 10
#' X <- matrix(mvrnorm(m*n,rep(0,m*n),diag(m*n)), nrow = n)
#' km_res <- gforce.kmeans(X,3)
#'
#' @useDynLib GFORCE kmeans_pp_R
gforce.kmeans <- function(X,K,R_only=FALSE){
  res <- NULL
  if(!R_only){
    X <- t(X)  # C implementation expects columns to be data points
    m <- nrow(X)
    n <- ncol(X)
    result <- .C(kmeans_pp_R,
                 X = as.double(X),
                 K = as.integer(K),
                 n = as.integer(n),
                 m = as.integer(m),
                 group_assignments = as.integer(rep(0,n)),
                 centers = numeric(K*m))
    res$clusters <- result$group_assignments
    res$centers <- matrix(result$centers,ncol=K)
  } else{
    R_res <- kmeanspp(X,K)
    res$clusters <- R_res$cluster
    res$centers <- R_res$centers
  }
  return(res)
}

# R implementation of K-means++ for test comparison
kmeanspp <- function(D, K) {
  d <- dim(D)[1]
  centers <- rep(0,K)
  distances <- matrix(rep(0,d *(K - 1)), ncol = K - 1)

  prob_dist <- rep(1, d)
  for (i in 1:(K - 1)) {
    centers[i] <- sample.int(d, 1, prob = prob_dist)
    distances[, i] <- colSums((t(D) - D[centers[i], ])^2)
    prob_dist <- distances[cbind(1:d, max.col(-distances[, 1:i, drop = FALSE]))]
  }
  centers[K] <- sample.int(d, 1, prob = prob_dist)
  res <- kmeans(D, D[centers, ])
  return(res)
}