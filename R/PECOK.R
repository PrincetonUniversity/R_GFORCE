
#' Estimates \eqn{\Gamma} for the PECOK SDP.
#' 
#' In particular, it returns in diagonal form the estimator
#' \eqn{\Gamma} used to construct the PECOK penalized covariance estimator.
#' @param X \eqn{n x d} matrix. Each row represents a realization of a
#' \eqn{d} dimensional random vector.
#' @param par logical expression. If \code{par == TRUE}, then a multi-threaded version
#' of the function is called. If \code{par == FALSE}, a single-threaded version is called.
#' @return The estimator \eqn{\Gamma} as a \eqn{d} dimensional numeric vector.
#'
#' @examples
#' dat <- gforce.generator(5,20,20,3)
#' gam_hat <- gforce.Gamma(dat$X)
#'
#' @useDynLib GFORCE v_measure
#' @useDynLib GFORCE v_measure_par
#' @export
gforce.Gamma <- function(X,par=FALSE){
  dims <- dim(X)
  n <- dims[1]
  d <- dims[2]
  
  gd <- rep(0,d)
  ips <- t(X)%*%X
  ips_diag <- matrix(diag(ips),ncol=1)
  ones_d <- matrix(rep(1,d),ncol=1)
  n_xc_xd <- (ones_d%*%t(ips_diag) + ips_diag%*%t(ones_d) - 2*ips)^0.5
  
  if(par){
    result <- .C("v_measure_par",
                 IPS=as.double(ips),
                 n_xc_xd=as.double(n_xc_xd),
                 dimension=as.integer(d),
                 vm=numeric(d^2))
  } else{
    result <- .C("v_measure",
                 IPS=as.double(ips),
                 n_xc_xd=as.double(n_xc_xd),
                 dimension=as.integer(d),
                 vm=numeric(d^2))
  }
  
  Vs_upper <- matrix(result$vm,ncol=d)
  Vs <- Vs_upper + t(Vs_upper)
  for(a in 1:d){
    # get neighbors
    v_a <- Vs[a,]
    vamax <- max(v_a)
    v_a[a] <- 2*vamax
    ne1 <- which.min(v_a)
    v_a[ne1] <- 2*vamax
    ne2 <- which.min(v_a)
    gd[a] <- (1/n)*(ips[a,a] + ips[ne1,ne2] - ips[a,ne1] - ips[a,ne2])
  }
  
  return(gd)
}
