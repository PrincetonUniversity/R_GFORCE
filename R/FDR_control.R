#' FDR Control Procedure.
#'
#' Performs the Banjamini and Yeuketli FDR control procedure. As input it takes a symmetric matrix of
#' test statistics with standard normal null distributions.
#' 
#' @param test_stats symmetric matrix of test statistics.
#' @param alpha alpha level for the FDR control procedure.
#' @export
gforce.FDR_control <- function(test_stats,alpha) {
    # set up
    K <- nrow(test_stats)
    num_hypotheses <- (K^2 - K)/2
    N_BY <- sum(1/(1:num_hypotheses))
    beta <- alpha / (2*num_hypotheses*N_BY)
    
    # create absolute values of test_stats for
    # procedure
    test_stats_step_up <- abs(test_stats)
    for(i in 1:K){
        for(j in 1:i) {
            test_stats_step_up[i,j] <- 0
        }
    }
    # get descending tau value possibilities
    tau_levels <- sort(test_stats_step_up,decreasing=TRUE)
    tau_levels <- tau_levels[1:num_hypotheses]
    #R_tau_hat <- num_hypotheses
    # found_min_tau <- 0
    # tau_hat <- tau_levels[R_tau_hat]
    # while(!found_min_tau && R_tau_hat > 1){
    #     R_tau_hat <- R_tau_hat - 1
    #     tau_hat <- tau_levels[R_tau_hat]
    #     p_tau <- beta * R_tau_hat
    #     q_tau <- abs(qnorm(p_tau))
    #     if(tau_hat > q_tau){
    #         found_min_tau <- 1
    #         # R_tau_hat <- R_tau_hat + 1
    #         # tau_hat <- tau_levels[R_tau_hat]
    #     }
    # }

    # reject greater than or equal to tau
    R_tau_hat <- 0
    tau_hat <- Inf
    q_tau_hat <- Inf
    for(R_tau_new in 1:num_hypotheses){
        min_tau_new <- tau_levels[R_tau_new]
        p_tau <- beta*R_tau_new
        q_tau <- abs(qnorm(p_tau))
        if(min_tau_new >= q_tau){
            tau_hat <- min_tau_new
            R_tau_hat <- R_tau_new
            q_tau_hat <- q_tau
        }
    }

    # true_discoveries <- test_stats_step_up > tau_hat
    true_discoveries <- test_stats_step_up >= tau_hat

    res <- NULL
    res$reject_null <- true_discoveries
    res$R_tau_hat <- R_tau_hat
    res$tau_hat <- tau_hat

    return(res)
}


#' Convert confidence intervals to equivalent test statistics.
#'
#' Can convert a 4D array encoding the confidence intervals for a precision
#' matrix to standard normal test-statistics.
#' 
#' @param conf_ints symmetric matrix of confidence intervals.
#' @param alpha alpha level of the confidence intervals.
#' @export
gforce.confint2test <- function(conf_ints,alpha) {
    K <- nrow(conf_ints)
    test_stats <- matrix(rep(0,K^2),nrow=K)

    z_alpha <- qnorm(1 - (alpha/2))

    for(i in 1:K){
        for(j in 1:K){
            v_ij <- conf_ints[i,j,3] - conf_ints[i,j,1]
            v_ij <- v_ij / (2*z_alpha)
            test_stats[i,j] <- conf_ints[i,j,2]/v_ij #/(conf_ints[i,i,2]*v_ij)
        }
    }

    return(test_stats)
}