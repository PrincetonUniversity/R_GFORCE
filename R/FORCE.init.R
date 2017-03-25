# FUNCTION TO CONSTRUCT INITIAL ITERATE AND SOLUTION FOR
# RENEGARS METHOD APPLIED TO KMEANS SDP K-MEANS SDP

#' Generate Random Initial Iterate and Values for FORCE Algorithm.
#' 
#' This function constructs a random full-rank solution to the Peng-Wei \eqn{K}-means SDP
#' and also gives an initial iterate \eqn{X0} for the FORCE algorithm.
#' @param D \eqn{d x d} numeric array. Specifies the objective function of the \eqn{K}-means SDP.
#' @param K number of groups.
#' @param s numeric value. Indicates the initial mixing of \code{opt_estimate} and a random full rank solution.
#' @param opt_estimate \eqn{d x d} numeric array. Specifies a guess of a feasible solution and is used to construct the initial feasible solution.
#' @return An object with following components
#' \describe{
#' \item{\code{E}}{Strictly feasible solution.}
#' \item{\code{X0}}{Strictly feasible initial iterate.}
#' \item{\code{E_val}}{Objective value of SDP at \eqn{E}.}
#' \item{\code{X0_val}}{Objective value of SDP at \eqn{X_0}.}
#' }
#'
#'
#' @examples
#' K <- 5
#' n <- 50 
#' d <- 50
#' m <- 3
#' s <- 0.25
#' dat <- gforce.generator(K,d,n,m,graph='DeltaC')
#' sh <- t(dat$X)%*%dat$X / n
#' gh <- gforce.Gamma(dat$X,par=TRUE)
#' D <- diag(gh)-sh
#' opt_guess <- kmeans(sh,K)$cluster
#' init_vals <- gforce.FORCE.init(D,K,s)
#'
#'
#' @useDynLib GFORCE FORCE_initialization_R
#' @export
gforce.FORCE.init <- function(D,K,s,opt_estimate,R_only=TRUE) {
  res <- NULL
  if(!R_only){
    d <- dim(D)[1]
    C_result <- .C("FORCE_initialization_R",
                 D=as.double(D),
                 s=as.double(s),
                 d=as.integer(d),
                 K=as.integer(K),
                 opt_estimate=as.double(opt_estimate),
                 E=numeric(d^2),
                 X0=numeric(d^2),
                 E_obj=as.double(0.0),
                 X0_obj=as.double(0.0))
    X0 <- C_result$X0
    E <- C_result$E
    dim(E) <- c(d,d)
    dim(X0) <- c(d,d)
    res$E <- E
    res$X0 <- X0
    res$E_val <- C_result$E_val
    res$X0_val <- C_result$X0_val
  } else{
    res <- renegar_start(D,K,s,opt_estimate)
  }
  return(res)
}

renegar_start <- function(diff,K,s,opt_estimate) {
  d <- dim(diff)[1]
  i <- 0
  full_rank_base <- gforce.full_rank_feasible(d,K)
  E <- matrix(rep(0,d^2),ncol=d)
  while(i < d){
    rotation <- sample(d)
    E <- E + full_rank_base[rotation,rotation]
    i <- i+1
  }

  E <- E / d
  E_val <- sum(diff*E)

  i <- 0
  M2 <- matrix(rep(0,d^2),ncol=d)
  while(i < d) {
    rotation <- sample(d)
    M2 <- M2 + full_rank_base[rotation,rotation]
    i <- i+1
  }

  M2 <- M2 / d
  M2_val <- sum(diff*M2)

  if(M2_val > E_val){
    X0 <- E
    X0_val <- E_val
    E <- M2
    E_val <- M2_val
  } else{
    X0 <- M2
    X0_val <- M2_val
  }
  res <- NULL
  res$E <- E
  res$E_val <- E_val
  res$X0 <- (1-s)*X0 + s*opt_estimate
  res$X0_val <- sum(diff*res$X0)
  return(res)
}