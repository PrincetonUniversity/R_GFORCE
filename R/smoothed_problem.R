# functions for evaluating the objective value and gradient
# of the smoothed problem from Renegar 2014

smoothed_gradient <- function(X,E,ESI,mu){
  grad <- NULL
  X_EVEV <- eigen(ESI%*%X%*%ESI,symmetric=TRUE)
  V_X <- Re(X_EVEV$vectors)
  X_eigs <- Re(X_EVEV$values)
  S_eigs <- X/E
  lambda_min <- min(min(X_eigs),min(S_eigs))
  X_eigs <- X_eigs - lambda_min
  d_X_new <- exp(-X_eigs / mu)
  S_eigs <- S_eigs - lambda_min
  GS <- exp(-S_eigs / mu)
  scale_factor <- sum(d_X_new) + sum(GS)
  grad$GX <- V_X %*% diag(d_X_new) %*% t(V_X) / scale_factor
  grad$GS <- GS / scale_factor
  return(grad)
}

smoothed_objective <- function(X,E,ESI,mu){
  X_EVEV <- eigen(ESI%*%X%*%ESI,symmetric=TRUE,only.values = TRUE)
  X_eigs <- Re(X_EVEV$values)
  S_eigs <- X/E
  lambda_min <- min(min(X_eigs),min(S_eigs))
  X_eigs <- X_eigs - lambda_min
  S_eigs <- S_eigs - lambda_min
  X_eig_sum <- sum(exp(-X_eigs/mu))
  S_eig_sum <- sum(exp(-S_eigs/mu))
  obj_val <- -mu*log(X_eig_sum + S_eig_sum) + lambda_min
  
  res <- NULL
  res$lambda_min <- lambda_min
  res$objective_value <-obj_val
  return(res)
}