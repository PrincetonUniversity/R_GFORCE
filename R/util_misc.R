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