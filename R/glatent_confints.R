
#' Confidence Intervals for Estimation in G-Latent Models.
#'
#' Estimate the precision matrix and construct confidence intervals for the latent and group averages graphs. The user
#' can either provide an estimate of the relevant covariance matrix or the data and cluster assignments. If cross validation
#' is selected, the data and clusters must be provided.
#' 
#' @param K number of clusters.
#' @export
gforce.glatent_confints <- function(C_hat=NULL,X_vals=NULL,clusters=NULL,alpha = 0.05,graph='latent',use_cv = FALSE,cv_opts = NULL,lambda1=NULL,lambda2=NULL,use_scio_package=TRUE) {
  res <- NULL
  # if user specified C_hat
  if(!is.null(C_hat)){
    if(graph == 'latent') {
      # cross validation cannot be used
      ;
    } else if(graph == 'averages') {
      # cross validation cannot be used
      ;
    } else {
      stop('gforce.glatent_confints -- You must specify either latent or averages graphs.')
    }

  } else if(!is.null(clusters) && !is.null(X_vals)) {
    if(graph == 'latent') {
      if(use_cv){
        res <- latent_confidence_intervals_all_cv(X_vals,clusters,alpha,use_scio_package,cv_opts)
      } else {
        ;
      }
    } else if(graph == 'averages') {
      if(use_cv){
        res <- averages_confidence_intervals_all_cv(X_vals,clusters,alpha,use_scio_package,cv_opts)
      } else {
        ;
      }
    } else {
      stop('gforce.glatent_confints -- You must specify either latent or averages graphs.')
    }


  } else {
    stop('gforce.glatent_confints -- You must specify one of C_hat or clusters and X.')
  }
  return(res)
}

#' Default Cross Validation Options for Confidence Intervals.
#'
#' 
#' @export
gforce.glatent_confints.cv_defaults <- function() {
  res <- NULL
  res$grid_density <- 5
  res$lambda1_min_exp <- -8
  res$lambda1_max_exp <- 0
  res$lambda2_min_exp <- -8
  res$lambda2_max_exp <- 0
  res$num_folds <- 5
  return(res)
}


latent_confidence_intervals_all_cv <- function(X_vals,group_assignments,alpha,use_scio_package,cv_opts=NULL) {
  # check for cv_opts
  if(is.null(cv_opts)){
    cv_opts <- gforce.glatent_confints.cv_defaults()
  }

  # ensure group labels ordered properly
  group_assignments <- order_group_assignments(group_assignments)
  
  # create C hat estimator
  n <- nrow(X_vals)
  groups <- unique(group_assignments)
  K <- length(groups)
  res <- array(0,c(K,K,3))
  sig_hat <- (t(X_vals)%*%X_vals)/n
  gam_hat <- glatent_Gamma_hat(sig_hat,group_assignments)
  Chat <- latent_spectrum_check(C_hat(sig_hat,gam_hat,group_assignments))
  
  #estimate theta *columnwise*
  theta_hat <- NULL
  if(!use_scio_package){
    theta_hat <- array(0,c(K,K))
    for(k in 1:K){
      objective_function <- function(ch,tha) scio_objective_function(ch,tha,k)
      function_solver <- function(ch,lambda) gforce.scio(ch,lambda,k)
      lambda1 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,cv_opts$lambda1_min_exp,cv_opts$lambda1_max_exp,cv_opts$num_folds)
      theta_hat[,k] <- gforce.scio(Chat,lambda1,k)
    }
    
    
  } else {
    lambda1 <- cv_lambda_selection_scio(X_vals,group_assignments,cv_opts$num_folds,1,20,0.7)
    theta_hat <- scio_package(Chat,lambda1)
  }

  #estimate v *row-wise*
  v_hat <- array(0,c(K,K))
  for(t in 1:K){
    objective_function <- function(ch,tha) max(abs(ch%*%tha))
    function_solver <- function(ch,lambda) v_hat(ch,t,lambda)
    lambda2 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,cv_opts$lambda2_min_exp,cv_opts$lambda2_max_exp,cv_opts$num_folds)
    v_hat[t,] <- v_hat(Chat,t,lambda2)
  }

  #estimate \tilde \theta_t,k *row-wise*
  for(t in 1:K){
    res[t,,2] <- decorrelated_estimator_t(Chat,theta_hat,v_hat[t,],t)
  }
  
  #construct upper and lower confidence bounds
  z_score <- qnorm(1 - (alpha/2))

  latent_variances <- latent_graph_variances(theta_hat,v_hat,gam_hat,group_assignments)

  for(k in 1:K){
    for(t in 1:K){
      conf_range <- sqrt(latent_variances[t,k])*z_score/sqrt(n)
      res[t,k,1] <- res[t,k,2] - conf_range
      res[t,k,3] <- res[t,k,2] + conf_range
    }
  }
  
  return(res)
}


# construct decorrelated confidence intervals for all entries
# returns 3 dimensional array. Last dimension is lower confidence bound, deco estimate, upper confidence bound
latent_confidence_intervals_all <- function(Chat,n,alpha,lambda1,lambda2) {
  K <- nrow(Chat)
  res <- array(0,c(K,K,3))
  
  #estimate theta *columnwise*
  theta_hat <- array(0,c(K,K))
  for(k in 1:K){
    theta_hat[,k] <- gforce.scio(Chat,lambda1,k)
  }
  
  #estimate v *row-wise*
  v_hat <- array(0,c(K,K))
  for(t in 1:K){
    v_hat[t,] <- v_hat(Chat,t,lambda2)
  }
  
  #estimate \tilde \theta_t,k *row-wise*
  for(t in 1:K){
    res[t,,2] <- decorrelated_estimator_t(Chat,theta_hat,v_hat[t,],t)
  }
  
  #construct upper and lower confidence bounds
  z_score <- qnorm(1 - (alpha/2))
  for(t in 1:K){
    for(k in 1:K){
      variance <- NULL
      if(t == k){
        c_tt <- Chat[t,t]
        t_kt <- theta_hat[k,t]
        t_kk <- theta_hat[k,k]
        t_tt <- theta_hat[t,t]
        variance <- c_tt^2 + t_kk*( 3*c_tt - 3*c_tt^2*t_tt + c_tt^3*t_tt^2 )
      }else{
        c_tt <- Chat[t,t]
        t_kk <- theta_hat[k,k]
        t_tt <- theta_hat[t,t]
        variance <- 1 +t_kk*(3*c_tt - 3*c_tt^2*t_tt + c_tt^3*t_tt^2)
      }
      conf_range <- sqrt(variance)*z_score/sqrt(n)
      res[t,k,1] <- res[t,k,2] - conf_range
      res[t,k,3] <- res[t,k,2] + conf_range
    }
  }
  
  return(res)
}


averages_confidence_intervals_all_cv <- function(X_vals,group_assignments,alpha,use_scio_package,cv_opts=NULL) {
  # check for cv_opts
  if(is.null(cv_opts)){
    cv_opts <- gforce.glatent_confints.cv_defaults()
  }

  # ensure group labels ordered properly
  group_assignments <- order_group_assignments(group_assignments)
  
  # create C hat estimator
  n <- nrow(X_vals)
  groups <- unique(group_assignments)
  K <- length(groups)
  res <- array(0,c(K,K,3))
  sig_hat <- (t(X_vals)%*%X_vals)/n
  s_hat <- S_hat(sig_hat,group_assignments)
  
  #estimate theta *columnwise*
  xi_hat <- array(0,c(K,K))
  for(k in 1:K){
    objective_function <- function(ch,tha) scio_objective_function(ch,tha,k)
    function_solver <- function(ch,lambda) gforce.scio(ch,lambda,k)
    lambda1 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,cv_opts$lambda1_min_exp,cv_opts$lambda1_max_exp,cv_opts$num_folds)
    xi_hat[,k] <- gforce.scio(s_hat,lambda1,k)
  }
  
  #estimate v *row-wise*
  v_hat <- array(0,c(K,K))
  for(t in 1:K){
    objective_function <- function(ch,tha) max(abs(ch%*%tha))
    function_solver <- function(ch,lambda) v_hat(ch,t,lambda)
    lambda2 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,cv_opts$lambda2_min_exp,cv_opts$lambda2_max_exp,cv_opts$num_folds)
    v_hat[t,] <- v_hat(s_hat,t,lambda2)
  }
  
  #estimate \tilde \theta_t,k *row-wise*
  for(t in 1:K){
    res[t,,2] <- decorrelated_estimator_t(s_hat,xi_hat,v_hat[t,],t)
  }
  
  #construct upper and lower confidence bounds
  z_score <- qnorm(1 - (alpha/2))
  for(t in 1:K){
    for(k in 1:K){
      variance <- (xi_hat[t,k]^2 + xi_hat[t,t]*xi_hat[k,k])#/(xi_hat[t,t]^2)
      conf_range <- sqrt(variance)*z_score/sqrt(n)
      res[t,k,1] <- res[t,k,2] - conf_range
      res[t,k,3] <- res[t,k,2] + conf_range
    }
  }
  
  return(res)
}

