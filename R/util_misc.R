# functions to build the constraint matrices for the projection
# onto the constraint set and for the dual construction

I_ab <- function(d,a,b) {
  Iab <- matrix(rep(0,d^2),ncol=d)
  Iab[a,b] <- 0.5
  Iab[b,a] <- 0.5
  return(Iab)
}

R_a <- function(d,a) {
  j <- matrix(rep(1,d),ncol=1)
  ea <- matrix(rep(0,d),ncol=1)
  ea[a] <- 1
  Ra <- ea%*%t(j) + j%*%t(ea)
  return(Ra)
}

dual_Q_feasible <- function(D,sol,y_a,y_T) {
  d <- ncol(D)
  Y_ab <- matrix(rep(0,d^2),ncol=d)
  for(i in 1:d){
    for(j in 1:i) {
      if(sol[i] != sol[j]){
        Y_ab[i,j] <- D[i,j] + y_a[i] + y_a[j]
        Y_ab[j,i] <- D[i,j] + y_a[i] + y_a[j]
      }
    }
  }

  Q <- D + y_T * diag(d)
  for(i in 1:d) {
    Q <- Q + y_a[i]*R_a(d,i)
  }
  Q <- Q - Y_ab

  res <- NULL
  res$Q <- Q
  res$Y_ab <- Y_ab
  return(res)
}

L_curve_criterion <- function(x) {
  res <- NULL

  d <- length(x)
  a <- c(1,x[1])
  b <- c(d,x[d])
  n <- b-a
  n <- n / sqrt(sum(n^2))

  dists <- rep(0,d)
  for(i in 1:d){
    c <- c(i,x[i])
    cproj_l <- sum((a-c)*n)
    cperp <- (a-c) - cproj_l*n
    dists[i] <- sqrt(sum(cperp^2))
  }
  res$dists <- dists
  res$max <- which.max(dists)
  return(res)
}

hclust_B <- function(BZ) {
  res <- NULL
  dBZ <- dist(BZ)
  hc <- hclust(dBZ)
  d <- ncol(BZ)

  MSEs <- rep(0,d)
  for(k in 2:d){
    cc <- cutree(hc,k=k)
    cc_mat <- gforce.clust2mat(cc)
    MSEs[k] <- 0.5*sum(cc_mat*as.matrix(dBZ))
  }

  lc <- L_curve_criterion(MSEs[2:(d-1)])
  res$K <- lc$max + 1
  res$clusters <- cutree(hc,k=res$K)
  res$MSE <- MSEs
  return(res)
}