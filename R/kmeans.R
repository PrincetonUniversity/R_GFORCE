# Interface to C implementation of K-means++ and Lloyd's Algorithm

# centers are the columns

#' @useDynLib GFORCE kmeans_pp_R
gforce.kmeans <- function(X,K){
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
  res <- NULL
  res$clusters <- result$group_assignments
  res$centers <- matrix(result$centers,ncol=K)
  return(res)
}