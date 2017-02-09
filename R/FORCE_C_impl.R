# Interface to C solver for K-means SDP

# FORCE solves the minimization form of the K-means SDP
# e.g. D can be a distances matrix or negated covariance
# INPUTS: Needs to have D and K. Optional force_opts
# optional initial feasible, optional E
# OUTPUTS: 

#' @useDynLib GFORCE primal_dual_adar_R
gforce.FORCE <- function(D,K,force_opts = NULL,D_Kmeans = NULL, X0 = NULL, E = NULL) {
    d <- ncol(D)

    if(is.null(force_opts)){
        force_opts <- pgd_adar_default_options_fixed(d,K)
    }

    if(is.null(D_Kmeans)) {
        D_Kmeans <- D
    }

    km_res <- kmeanspp(-D_Kmeans,K)
    km_sol <- B_hat(km_res)

    if(is.null(X0) && is.null(E)){
        ren_start_res <- renegar_start(D,K,force_opts$initial_mixing,km_sol)
        X0 <- ren_start_res$X1
        E <- ren_start_res$E 
    } else if(is.null(X0)){
        stop('FORCE -- Either specify both X0 and E or neither\r\n')
    } else if(is.null(E)){
        stop('FORCE -- Either specify both X0 and E or neither\r\n')
    }

    E_EVEV <- eigen(E)
    E_V <- E_EVEV$vectors
    E_D <- diag(E_EVEV$values)
    E_sqrt <- E_V%*%(E_D^(0.5))%*%t(E_V)
    ESI <- solve(E_sqrt)

    C_result <- .C(primal_dual_adar_R,
            D = as.double(D),
            D_Kmeans = as.double(D_Kmeans),
            E = as.double(E),
            ESI = as.double(ESI),
            X0 = as.double(X0),
            d = as.integer(d),
            K = as.integer(K),
            verbosity = as.integer(force_opts$verbose),
            kmeans_iter = as.integer(force_opts$kmeans_iter),
            dual_frequency = as.integer(force_opts$dual_frequency),
            max_iter = as.integer(force_opts$max_iter),
            finish_pgd = as.integer(force_opts$finish_pgd),
            number_restarts = as.integer(length(force_opts$restarts)),
            restarts = as.integer(force_opts$restarts),
            alpha = as.double(force_opts$alpha),
            eps_obj = as.double(force_opts$eps_obj),
            Z_T = numeric(d^2),
            B_Z_T = numeric(d^2),
            Z_T_lmin = as.double(1.0),
            Z_best = numeric(d^2),
            B_Z_best = numeric(d^2),
            Z_best_lmin = as.double(1.0),
            B_Z_T_opt_val = as.double(1.0),
            B_Z_best_opt_val = as.double(1.0),
            km_opt_val = as.double(1.0),
            km_best = as.integer(1:d),
            km_best_time = as.integer(1),
            km_iter_best = as.integer(1),
            km_iter_total = as.integer(1),
            dc = as.integer(1),
            dc_time = as.integer(1),
            dc_grad_iter = as.integer(1),
            grad_iter_best = as.integer(1),
            grad_iter_best_time = as.integer(1),
            total_time = as.integer(1))

    # Build Result
    res <- NULL
    res$X0 <- X0
    res$E <- E
    res$Z_T <- matrix(C_result$Z_T,ncol=d)
    res$B_Z_T <- matrix(C_result$B_Z_T,ncol=d)
    res$B_Z_T_opt_val <- C_result$B_Z_T_opt_val
    res$Z_best <- matrix(C_result$Z_best,ncol=d)
    res$B_Z_best <- matrix(C_result$B_Z_best,ncol=d)
    res$B_Z_best_opt_val <- C_result$B_Z_best_opt_val
    res$km_best <- C_result$km_best
    res$km_opt_val <- C_result$km_opt_val
    res$B_km <- B_hat(res$km_best)
    res$km_best_time <- C_result$km_best_time
    res$km_iter_best <- C_result$km_iter_best
    res$km_iter_total <- C_result$km_iter_total
    res$dual_certified <- C_result$dc
    res$dual_certified_grad_iter <- C_result$dc_grad_iter
    res$grad_iter_best <- C_result$grad_iter_best
    res$grad_iter_best_time <- C_result$grad_iter_best_time
    res$total_time <- C_result$total_time

    return(res)
}


# Solve the PECOK SDP using the FORCE algorithm
# Same outputs as FORCE, but also will include the
# matrix D
gforce.PECOK <- function(K, X=NULL, D=NULL, force_opts = NULL) {

}