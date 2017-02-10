# FUNCTIONS FOR CLUSTERING AND METRICS

check_perfect_recovery <- function(original,recovered){
  group_ids_orig <- unique(original)
  group_ids_recov <- unique(recovered)
  K1 <- length(group_ids_orig)
  K2 <- length(group_ids_recov)
  if(K1 != K2){
    return(FALSE)
  }
  K <- K1
  original_groups <- list()
  recovered_groups <- list()
  for(i in 1:K){
    original_groups[[i]]  <- which(original == group_ids_orig[i])
    recovered_groups[[i]]  <- which(recovered == group_ids_recov[i])
  }
  same <- TRUE
  for(i in 1:K){
    orig_group <- original_groups[[i]]
    found_match <- FALSE
    for(j in 1:K){
      # check for match in recovered group j
      recov_group <- recovered_groups[[j]]
      if(length(orig_group) == length(recov_group)){
        if(length(orig_group) == sum(orig_group %in% recov_group)){
          found_match <- TRUE
        }
      }
    }
    if(!found_match){
      same <- FALSE
    }
  }
  
  return(same)
}

purity_measure <- function(ga_hat,ga){
  group_ids_ga <- unique(ga)
  group_ids_ga_hat <- unique(ga_hat)
  K <- length(group_ids_ga_hat)
  K2 <- length(unique(group_ids_ga))
  n <- length(ga_hat)
  p <- 0
  for(k in 1:K){
    m <- 0
    for(kp in 1:K2){
      group_ga_hat_k <- which(ga_hat == group_ids_ga_hat[k])
      group_ga_kp <- which(ga == group_ids_ga[kp])
      int_size <- sum(group_ga_hat_k %in% group_ga_kp)
      m <- max(m,int_size)
    }
    p <- p + m
  }
  p <- p/n
  return(p)
}

# R implementation of K-means++ for test comparison
kmeanspp <- function(D, K) {
  d <- dim(D)[1]
  centers <- rep(0,K)
  distances <- matrix(rep(0,d *(K - 1)), ncol = K - 1)

  prob_dist <- rep(1, d)
  for (i in 1:(K - 1)) {
    centers[i] <- sample.int(d, 1, prob = prob_dist)
    distances[, i] <- colSums((t(D) - D[centers[i], ])^2)
    prob_dist <- distances[cbind(1:d, max.col(-distances[, 1:i, drop = FALSE]))]
  }
  centers[K] <- sample.int(d, 1, prob = prob_dist)
  res <- kmeans(D, D[centers, ])
  return(res$cluster)
}

misclassified_points <- function(ga_hat,ga){
  group_ids_ga <- unique(ga)
  group_ids_ga_hat <- unique(ga_hat)
  K <- length(group_ids_ga_hat)
  K2 <- length(unique(group_ids_ga))
  n <- length(ga_hat)
  p <- 0
  misclassified_points <- list()
  for(k in 1:K){
    m <- 0
    misclassified_points_from_grp <- c()
    group_ga_hat_k <- which(ga_hat == group_ids_ga_hat[k])
    for(kp in 1:K2){
      group_ga_kp <- which(ga == group_ids_ga[kp])
      corr_class <- group_ga_hat_k %in% group_ga_kp
      int_size <- sum(corr_class)
      if(int_size > m){
        misclassified_points_from_grp <- setdiff(group_ga_hat_k,group_ga_kp[corr_class])
        m <- int_size
      }
    }
    misclassified_points[[k]] <- misclassified_points_from_grp
  }
  return(misclassified_points)
}

kmeans_repeater <- function(sig,num_repeat,ga){
  source(paste(DIR_CONVEX_KMEANS,'clustering.R',sep='/'))
  av_purity <- 0
  percent_perfect <- 0
  K <- length(unique(ga))
  for(i in 1:num_repeat){
    cur_purity <- purity_measure(kmeanspp(sig,K),ga)
    av_purity <- av_purity + cur_purity
    if(cur_purity == 1){
      percent_perfect <- percent_perfect+1
    }
  }
  res <- NULL
  res$percent_perfect <- percent_perfect/num_repeat
  res$average <- av_purity/num_repeat
  return(res)
}

