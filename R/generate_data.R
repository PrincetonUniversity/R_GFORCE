# GENERATE RANDOM DATA ACCORDING TO A G-BLOCK LATENT VARIABLE MODEL

#' Data generator.
#' 
#' Generates \eqn{n} random samples from a \eqn{G}-Latent Variable Model. The
#' caller can specify the graph structure on the latent variables via
#' several parameters. The magnitude of the non-zero entries in the population
#' precision matrix can also be specified. Observed variables are assigned
#' uniformly at random to \eqn{K} groups with minimum size \eqn{m}.
#'
#' @param K number of clusters.
#' @param n number of samples.
#' @param d dimension of the observed random vector.
#' @param m minimal group size.
#' @param graph latent graph structure. Can be 'scalefree', 'hub', 'band' or 'DeltaC'.
#' @param num_hubs number of hubs in the latent graph. Ignored unless \code{graph == 'hub'}.
#' @param band_size size of bands in the latent graph. Ignored unless \code{graph=='band'}.
#' @param cov_gap_mult scales the size of \eqn{\delta C}. Ignored unless \code{graph == 'DeltaC'}.
#' @param error_base minimum variance of \eqn{E_i}.
#' @param error_add size of range of possible variances for \eqn{E_i}.
#' @param corr_value size of off diagonal entries in latent precision matrix.
#' @param normalize logical. If \code{normalize == TRUE}, the covariance matrix for the latent graph
#' will be normalized so that it is also a correlation matrix.
#' 
#' @export
gforce.generator <- function(K,d,n,m, graph = 'DeltaC', num_hubs=NULL, band_size = 3, cov_gap_mult=1.0, 
                            error_base = 0.25, error_add = 0.25, corr_value = 0.3, normalize = TRUE) {
  require(MASS)
  res <- NULL

  # generate error variances
  gamma_star <- runif(d,max=error_add)
  res$gamma_star <- gamma_star + error_base

  # build latent covariance structure
  if(graph == 'DeltaC'){
    gamma_inf <- max(res$gamma_star)
    inter <- sqrt(log(d)/(m*n)) + sqrt(d / (n*m^2)) + log(d)/n + d/(m*n)
    zeta <- cov_gap_mult*gamma_inf*inter
    res$Cstar <- random_covariance(K,zeta)
    res$Theta_star <- solve(res$Cstar)

  } else if(graph == 'scalefree') {
    A <- matrix(rep(0,K^2), nrow=K)
    A[1,2] <- 0.3
    A[2,1] <- 0.3
    for(i in 3:K){
      # get current probability distribution
      probd <- colSums(A)
      probd <- probd / sum(probd)
      
      #add new edge and node
      ne <- sample(x = 1:K, 1, replace = T, prob = probd)
      A[i,ne] <- 0.3
      A[ne,i] <- 0.3
    }
    res$Theta_star <- A + (abs(min(eigen(A)$values)) + 0.2)*diag(K)

  } else if(graph == 'band') {
    A <- matrix(rep(0,K^2), nrow=K)
    for(i in 1:band_size){
      diag(A[1:(K-i),(1+i):K]) <- corr_value
      diag(A[(1+i):K,1:(K-1)]) <- corr_value
    }
    res$Theta_star <- A + (abs(min(eigen(A)$values)) + 0.2)*diag(K)

  } else if(graph == 'hub') {
    if(K %% num_hubs != 0){
      stop('Number of Hubs must divide K')
    }
    A <- matrix(rep(0,K^2), nrow=K)
    group_size <- K / num_hubs
    for(i in 1:num_hubs){
      cur_hub <- ((i-1)*group_size + 1):(i*group_size)
      hub_center <- cur_hub[1]
      A[cur_hub[-1],hub_center] = corr_value
      A[hub_center,cur_hub[-1]] = corr_value
    }
    res$Theta_star <- A + (abs(min(eigen(A)$values)) + 0.2)*diag(K)

  } else {
    stop(sprintf('gforce.generator -- Graph type %c not supported.',graph))
  }

  if(graph != 'DeltaC'){
    if(normalize){
      res$Cstar <- cov2cor(solve(res$Theta_star))
      res$Theta_star <- solve(res$Cstar)
    } else {
      res$Cstar <- solve(res$Theta_star)
    }
  }

  #group assignments
  res$group_assignments <- generate_random_partition(K,d,m)

  # Generate Data
  res$E <- mvrnorm(n,rep(0,d),diag(res$gamma_star))
  res$Z <- mvrnorm(n,rep(0,K),res$Cstar)
  res$X <- matrix(rep(0,d*n),nrow=n)
  for(i in 1:K) {
    group_idx = which(res$group_assignments == i)
    res$X[,group_idx] = res$E[,group_idx] + res$Z[,i]
  }
  
  return(res)
}

#compute the covariance gap
delta_c <- function(CS) {
  K <- dim(CS)[1]
  dc = 1000000000000;
  for(i in 2:K){
    for(j in 1:(i-1)){
      difc <- CS[i,i] + CS[j,j] - 2*CS[i,j]
      dc <- min(difc,dc)
    }
  }
  return(dc)
}

random_covariance <- function(K,zeta) {
  min_eig <- 0
  dc <- 0
  tries <- 0
  C <- NULL
  while(((min_eig <= 0) || (dc < zeta)) && (tries < 20)) {
    a <- matrix(rnorm(K^2),ncol=K)
    C_new <- a%*%t(a)
    min_eig_new <- min(eigen(C_new)$values)
    dc_new <- delta_c(C_new)
    if(dc_new > dc && min_eig_new > 0) {
      C <- C_new
      dc <- dc_new
    }
    tries <- tries + 1
  }
  C <- (zeta/dc)*C
  return(C)
}


generate_random_partition <- function(K,d,m){
  group_sum <- 0
  pre_alloc <- K*m
  group_sizes <- NULL
  while(group_sum != (d-pre_alloc)){
    group_sizes <- runif(K)
    group_sum = sum(group_sizes)
    group_sizes <- ((d - pre_alloc)/group_sum)*group_sizes
    group_sizes <- round(group_sizes)
    group_sum <- sum(group_sizes)
  }
  group_sizes <- group_sizes + m
  group_assignments <- rep(0,d)
  group_number <- 1
  group_counter <- 0
  for(i in 1:d) {
    if(group_counter < group_sizes[group_number]) {
      group_counter <- group_counter+1
      group_assignments[i] <- group_number
    }
    else{
      group_number <- group_number + 1
      group_counter <- 1
      group_assignments[i] <- group_number
    }
  }
  return(group_assignments)
}