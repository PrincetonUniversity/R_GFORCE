# FUNCTIONS TO CONSTRUCT OPTIMALITY CERTIFICATES FOR THE 
# K-MEANS SDP













# This function is used now for testing the C implementation
# construct dual solutions to kmeans SDP relaxation
# ga_hat -- proposed primal integer solution
# D -- objective value function in **maximization** version of the relaxation
dual_solution <- function(ga_hat,D,eps = 0.01,eps2 = 10^-7,Y_T_min = 0.01){

  # initialization and pre-processing
  # source(paste(DIR_CONVEX_KMEANS,'constraint_matrices.R',sep='/'))
  # source(paste(DIR_CONVEX_KMEANS,'B_hat.R',sep='/'))
  return_dual <- NULL
  group_ids <- unique(ga_hat)
  K <- length(group_ids)
  d <- dim(D)[1]
  group_idxs <- list()
  group_sizes <- rep(0,K)
  group_sums <- rep(0,K)
  for(i in 1:K){
    ga_hat_k <- which(ga_hat == group_ids[i])
    group_idxs[[i]] <- ga_hat_k
    group_sizes[i] <- length(ga_hat_k)
    group_sums[i] <- sum(D[ga_hat_k,ga_hat_k])
  }
  
  # get primal optimal
  B_opt <- B_hat(ga_hat)
  primal_value <- sum(D*B_opt)
  
  # construct dual optimal solution by binary search
  # on Y_T
  M <- matrix(rep(0,d*(d+1)),ncol=d)
  b_base <- matrix(rep(0,d+1),ncol=1)
  
  for(k in 1:K){
    g_k_idx <- group_idxs[[k]]
    g_k_size <- group_sizes[k]
    M[g_k_idx,g_k_idx] <- matrix(rep(1,g_k_size^2),ncol=g_k_size) + g_k_size*diag(g_k_size)
    if(g_k_size == 1){
      b_base[g_k_idx] <- D[g_k_idx,g_k_idx]
    } else{
      b_base[g_k_idx] <- colSums(D[g_k_idx,g_k_idx])      
    }
  }
  M[d+1,1:d] <- 2*rep(1,d)
  b_base[d+1] <- primal_value
  
  Y_T_max <- abs(primal_value)
  dual_best <- NULL
  found_feasible <- FALSE
  
  # perform the search
  while(Y_T_max / Y_T_min > 1+eps){
    # construct b vector for this iteration
    Y_T_new <- (Y_T_max + Y_T_min) / 2
    b_new <- rep(0,d+1)
    b_new[1:d] <- b_base[1:d] - Y_T_new
    b_new[d+1] <- b_base[d+1] - Y_T_new*K
    
    # solve linear system
    Y_a_new <- qr.solve(M,b_new)
    
    # update current solution
    A <- Y_T_new*diag(d)
    for(a in 1:d){
      A <- A + Y_a_new[a]*R_a(d,a)
    }
    Y_ab_new <- A - D
    Y_ab_new[B_opt > 0] <- 0
    R_new <- A - D - Y_ab_new
    R_evals <- eigen(R_new,only.values=TRUE)$values
    em_new <- min(Re(R_evals))
    
    # update bounds on search
    if(em_new > -eps2){
      found_feasible <- TRUE
      dual_best$Y_T <- Y_T_new
      dual_best$Y_a <- Y_a_new
      dual_best$Y_ab <- Y_ab_new
      dual_best$R <- R_new
      Y_T_max <- Y_T_new
    } else {
      Y_T_min <- Y_T_new
    }
  }
  
  # clean up final solution if needed to ensure feasibility
  if(found_feasible){
    return_dual$Y_a <- dual_best$Y_a
    return_dual$Y_ab <- dual_best$Y_ab
    return_dual$Y_ab[return_dual$Y_ab < 0] = 0
    A <- dual_best$Y_T*diag(d)
    for(a in 1:d){
      A <- A + return_dual$Y_a[a]*R_a(d,a)
    }
    R <- A - D - return_dual$Y_ab
    R_evals <- eigen(R,only.values=TRUE)$values
    return_dual$Y_T <- dual_best$Y_T - min(0,min(Re(R_evals)))
    return_dual$R <- R + (return_dual$Y_T - dual_best$Y_T)*diag(d)
    return_dual$primal_value <- primal_value
    return_dual$dual_value <- 2*sum(return_dual$Y_a) + K*return_dual$Y_T
  }
  
  return(return_dual)
}