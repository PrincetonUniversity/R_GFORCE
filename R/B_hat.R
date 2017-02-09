# construct solution to kmeans relaxation from group assignment vector

B_hat <- function(ga) {
  d <- length(ga)
  K <- length(unique(ga))
  Bhat <- matrix(rep(0,d^2),ncol=d)
  
  for(k in 1:K){
    idx <- which(ga == k)
    group_size <- length(idx)
    for(a in idx){
      for(b in idx){
        Bhat[a,b] <- 1/group_size
      }
    }
  }
  return(Bhat)
}