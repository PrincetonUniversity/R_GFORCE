# Hierarchical Clustering with cluster estimation.

#' Hierarchical Clustering with Estimation of \eqn{K}.
#'
#' Clusters input and estimates \eqn{K}.
#' @param X \eqn{n x m} matrix. Each row is treated as a point in \eqn{R^m}.
#' @param K integer. The number of clusters to group the data into.
#' @param R_only logical expression. If \code{R_only == FALSE}, then the included
#' native code implementation will be used. Otherwise, an R implementation is used.
#' @return Returns an object with the components:
#' \describe{
#' \item{\code{K}}{an estimate of the number of clusters.}
#' \item{\code{clusters}}{a \eqn{n} dimensional integer vector. Entry \eqn{i} to the cluster assignment of the data point given by row \eqn{i} of \code{X}.}
#' \item{\code{MSE}}{a \eqn{n} dimensional vector of the mean squared errors of each choice of \eqn{K}.}
#' }
#'
#' @examples
#' m <- 10 
#' n <- 10
#' X <- matrix(mvrnorm(m*n,rep(0,m*n),diag(m*n)), nrow = n)
#' hc_res <- gforce.hclust(X)
#'
#' @export
gforce.hclust <- function(X) {
  res <- NULL
  dX <- dist(X)
  hc <- hclust(dX)
  d <- ncol(X)

  MSEs <- rep(0,d)
  for(k in 2:d){
    cc <- cutree(hc,k=k)
    cc_mat <- gforce.clust2mat(cc)
    MSEs[k] <- 0.5*sum(cc_mat*as.matrix(dX))
  }

  lc <- L_curve_criterion(MSEs[2:(d-1)])
  res$K <- lc$max + 1
  res$clusters <- cutree(hc,k=res$K)
  res$MSE <- MSEs
  return(res)
}