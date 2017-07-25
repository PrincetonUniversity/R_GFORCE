# estimate a columnwise SCIO estimator like \hat\Theta_{\cdot K}
# see Liu and Luo (2012)

scio_estimator <- function(Chat, k, lambda, eps = 10^-18, max_iter = 1000) {
  # soft thresholding
  soft_threshold <- function(x){
    #find soft threshold
    st_x <- sign(x)*(abs(x)-lambda)
    return(st_x)
  }
  
  K <- nrow(Chat)
  scio_t <- NULL
  val_t <- 1
  val_tp1 <- 0
  t <- 0
  scio_tp1 <- rnorm(K) #rep(0,K)
  p <- 1
  while(abs(val_t-val_tp1) > eps && t < max_iter) {
    val_t <- val_tp1
    scio_t <- scio_tp1
    t <- t+1
    x_tp1_p <- scio_t[-p]%*%Chat[-p,p]
    s_tp1_p <- NULL
    if(p == k){
      s_tp1_p <- soft_threshold(1-x_tp1_p) / Chat[p,p]
    }else{
      s_tp1_p <- soft_threshold(-1*x_tp1_p) / Chat[p,p]
    }
    scio_tp1 <- scio_t
    scio_tp1[p] <- s_tp1_p
    if(p == K) {
      p <- 1
    } else{
      p <- p + 1
    }
    val_tp1 <- scio_objective_function(Chat,scio_tp1,k) + lambda*sum(abs(scio_tp1))
  }
  return(scio_tp1)
}

scio_objective_function <- function(Chat,theta,k){
  return((1/2)*(theta%*%Chat%*%theta) - theta[k])
}



scio_package <- function(Chat, lambda, eps=10^-6, max_iter=10000, sym=FALSE) {
  library(scio)
  scio_res <- scio(Chat,lambda,thr=eps, maxit = max_iter, sym=sym)
  return(scio_res$w)
}


scio_likelihood <- function(C, Cinv) {
  ot <- as.numeric(unlist(determinant(Cinv)))
  if (ot[2]<=0) warning("Precision matrix estimate is not positive definite!")
  tmp <- (sum(diag(C%*%Cinv))  - ot[1] - dim(Cinv)[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    return(Inf)
  }
}