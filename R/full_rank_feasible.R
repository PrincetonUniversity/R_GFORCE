# FUNCTION TO CONSTRUCT FULL RANK FEASIBLE SOLUTIONS TO 
# K-MEANS SDP

#' Strictly Feasible Solution to \eqn{K}-means SDP.
#' 
#' This function constructs a  full-rank solution to the Peng-Wei \eqn{K}-means SDP.
#' @param d dimension of the solution.
#' @param K number of groups.
#' @param R_only logical expression. If \code{R_only == TRUE}, then no native code is run.
#' @return A strictly feasible solution \eqn{E}.
#'
#'
#' @examples
#' K <- 5
#' n <- 50 
#' d <- 50
#' E <- gforce.full_rank_feasible(d,K)
#'
#' @useDynLib GFORCE full_rank_feasible_R
#' @export
gforce.full_rank_feasible <- function(d,K,R_only=FALSE) {
  if(!R_only){
    result <- .C("full_rank_feasible_R",
                 d=as.integer(d),
                 K=as.integer(K),
                 E=numeric(d^2))
    E <- result$E
    dim(E) <- c(d,d)
    return(E)
  } else {
    c <- floor(d/(K-1))
    E <- matrix(rep(0,d^2),ncol=d)
    for(i in 1:c){
      for(j in 1:(K-1)){
        p <- (i-1)*(K-1) + j
        E[p,p] <- E[p,p] + 1
      }
      ind_groups <- ((K-1)*(i-1) + 1):((K-1)*i)
      idx <- 1:d
      idx <- setdiff(idx,ind_groups)
      for(a in 1:(d-K+1)){
        for(b in 1:(d-K+1)){
          E[idx[a],idx[b]]<- E[idx[a],idx[b]] + 1/(d-K+1)
        }
      }
    }

    if((c - d/(K-1)) < 0){
      for(j in 0:(K-2)){
        p <- d-j
        E[p,p] <- E[p,p] + 1
      }
      ind_groups <- (d-K+2):d
      idx <- 1:d
      idx <- setdiff(idx,ind_groups)
      for(a in 1:(d-K+1)){
        for(b in 1:(d-K+1)){
          E[idx[a],idx[b]] <- E[idx[a],idx[b]] + 1/(d-K+1)
        }
      }
    }
    E <- E / ceiling(d/(K-1))
    return(E)
  }
  
}