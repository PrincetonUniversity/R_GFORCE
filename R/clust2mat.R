# Convert K-means solution to SDP solution

#' Convert K-means solution to partnership matrix.
#' 
#' @param clusters length \eqn{d} vector. Assigns each variable or data point to a cluster. Cluster names
#' can be numbers or strings.
#' @export
gforce.clust2mat <- function(clusters) {
  d <- length(clusters)
  group_names <- unique(clusters)
  K <- length(group_names)
  Bhat <- matrix(rep(0,d^2),ncol=d)
  
  for(k in 1:K){
    idx <- which(clusters == group_names[k])
    group_size <- length(idx)
    for(a in idx){
      for(b in idx){
        Bhat[a,b] <- 1/group_size
      }
    }
  }
  return(Bhat)
}