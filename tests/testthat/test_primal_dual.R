context('Test Primal Dual Algorithm')

#' @useDynLib GFORCE test_smoothed_gradient_S_base
test_that("GS_t Base",{

    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    initial_mixing <- 2/d
    km_res <- gforce.kmeans(-sh,K,R_only=TRUE)
    km_res <- km_res$clusters

    ren_start_res <- renegar_start(-sh,K,initial_mixing,km_res)
    X <- ren_start_res$X1
    E <- ren_start_res$E
    S_min_r <- 0

    result <- .C(test_smoothed_gradient_S_base,
                 X= as.double(X),
                 E = as.double(E),
                 GS_t = numeric(d^2),
                 d = as.integer(d),
                 S_min_r = as.double(S_min_r))
    GS_t <- matrix(result$GS_t,ncol=d)
    comp_GS_t <- X / E
    comp_Smin <- min(min(comp_GS_t))

    expect_equal(GS_t,comp_GS_t)
    expect_equal(comp_Smin,result$S_min_r)
    })

#' @useDynLib GFORCE test_smoothed_gradient_X_base
test_that("GX_t Base",{
    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    gh <- gforce.Gamma(dat$X)
    diff <- diag(gh) - sh
    initial_mixing <- 2/d
    km_res <- gforce.kmeans(-sh,K,R_only=TRUE)
    km_res <- km_res$clusters
    km_sol <- gforce.clust2mat(km_res)

    ren_start_res <- renegar_start(diff,K,initial_mixing,km_sol)
    X <- ren_start_res$X1
    E <- ren_start_res$E

    E_EVEV <- eigen(E)
    V <- E_EVEV$vectors
    D <- diag(E_EVEV$values)
    E_sqrt <- V%*%(D^(0.5))%*%t(V)
    ESI <- solve(E_sqrt)

    result <- .C(test_smoothed_gradient_X_base,
                 X= as.double(X),
                 ESI = as.double(ESI),
                 GX_t = numeric(d^2),
                 d = as.integer(d),
                 K = as.integer(K),
                 eig_vals = numeric(d))
    GX_t <- matrix(result$GX_t,ncol=d)
    E_X_E <- ESI%*%X%*%ESI
    X_EVEV <- eigen(E_X_E)
    comp_Xmin <- min(Re(X_EVEV$values))
    comp_E_X_E <- GX_t %*% diag(result$eig_vals) %*% t(GX_t)

    expect_equal(comp_Xmin,min(result$eig_vals))
    expect_equal(comp_E_X_E,E_X_E)
    })

#' @useDynLib GFORCE test_smoothed_gradient
test_that("Smoothed Gradient (GX_t, GS_t)",{
    set.seed(12345)
    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    gh <- gforce.Gamma(dat$X)
    diff <- diag(gh) - sh
    initial_mixing <- 2/d
    km_res <- gforce.kmeans(-sh,K,R_only=TRUE)
    km_res <- km_res$clusters
    km_sol <- gforce.clust2mat(km_res)

    ren_start_res <- renegar_start(diff,K,initial_mixing,km_sol)
    X <- ren_start_res$X1
    E <- ren_start_res$E

    E_EVEV <- eigen(E)
    V <- E_EVEV$vectors
    D <- diag(E_EVEV$values)
    E_sqrt <- V%*%(D^(0.5))%*%t(V)
    ESI <- solve(E_sqrt)

    mu <- 0.5*0.01/log(d)
    gradSX <- smoothed_gradient(X,E,ESI,mu)

    result <- .C(test_smoothed_gradient,
                 X= as.double(X),
                 E = as.double(E),
                 ESI = as.double(ESI),
                 d = as.integer(d),
                 K = as.integer(K),
                 mu = as.double(mu),
                 GX_t = numeric(d^2),
                 GS_t = numeric(d^2))
    GX_t <- matrix(result$GX_t,ncol=d)
    GS_t <- matrix(result$GS_t,ncol=d)

    expect_equal(GX_t,gradSX$GX)
    expect_equal(GS_t,gradSX$GS)
    })

#' @useDynLib GFORCE test_project_C_perpendicular
test_that("Projection onto C Perpendicular",{
    set.seed(12345)
    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    gh <- gforce.Gamma(dat$X)
    diff <- diag(gh) - sh
    initial_mixing <- 2/d
    km_res <- gforce.kmeans(-sh,K,R_only=TRUE)
    km_res <- km_res$clusters
    km_sol <- gforce.clust2mat(km_res)

    ren_start_res <- renegar_start(diff,K,initial_mixing,km_sol)
    X <- ren_start_res$X1
    E <- ren_start_res$E

    E_EVEV <- eigen(E)
    V <- E_EVEV$vectors
    D <- diag(E_EVEV$values)
    E_sqrt <- V%*%(D^(0.5))%*%t(V)
    ESI <- solve(E_sqrt)

    mu <- 0.5*0.01/log(d)
    gradSX <- smoothed_gradient(X,E,ESI,mu)

    result <- .C(test_project_C_perpendicular,
                 D= as.double(diff),
                 d = as.integer(d),
                 K = as.integer(K),
                 GX_t = as.double(gradSX$GX),
                 GS_t = as.double(gradSX$GS))
    Z_proj <- matrix(result$GX_t,ncol=d)

    comp_Z <- project_C_perpendicular(gradSX$GX, gradSX$GS,diff)

    expect_equal(comp_Z$Z_proj,Z_proj)
    })

#' @useDynLib GFORCE test_project_E
test_that("Projection Onto PSD Cone Border",{
    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    gh <- gforce.Gamma(dat$X)
    diff <- diag(gh) - sh
    initial_mixing <- 2/d
    km_res <- gforce.kmeans(-sh,K,R_only=TRUE)
    km_res <- km_res$clusters
    km_sol <- gforce.clust2mat(km_res)

    ren_start_res <- renegar_start(diff,K,initial_mixing,km_sol)
    X <- ren_start_res$X1
    E <- ren_start_res$E

    E_EVEV <- eigen(E)
    V <- E_EVEV$vectors
    D <- diag(E_EVEV$values)
    E_sqrt <- V%*%(D^(0.5))%*%t(V)
    ESI <- solve(E_sqrt)
    mu <- 0.5*0.01/log(d)
    res <- smoothed_objective(X,E,ESI,mu)
    lmin <- res$lambda_min
    result <- .C(test_project_E,
                 E= as.double(E),
                 X = as.double(X),
                 d = as.integer(d),
                 lmin = as.double(lmin),
                 Z = numeric(d^2))
    Z_proj <- matrix(result$Z,ncol=d)
    Z_proj_comp <- E + (1/(1-lmin))*(X - E)

    expect_equal(Z_proj_comp,Z_proj)
    })

#' @useDynLib GFORCE test_smoothed_objective
test_that("Smoothed Objective",{
    set.seed(12345)
    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    gh <- gforce.Gamma(dat$X)
    diff <- diag(gh) - sh
    initial_mixing <- 2/d
    km_res <- gforce.kmeans(-sh,K,R_only=TRUE)
    km_res <- km_res$clusters
    km_sol <- gforce.clust2mat(km_res)

    ren_start_res <- renegar_start(diff,K,initial_mixing,km_sol)
    X <- ren_start_res$X1
    E <- ren_start_res$E

    E_EVEV <- eigen(E)
    V <- E_EVEV$vectors
    D <- diag(E_EVEV$values)
    E_sqrt <- V%*%(D^(0.5))%*%t(V)
    ESI <- solve(E_sqrt)

    mu <- 0.5*0.01/log(d)
    comp_result <- smoothed_objective(X,E,ESI,mu)


    result <- .C(test_smoothed_objective,
                 X= as.double(X),
                 E = as.double(E),
                 ESI = as.double(ESI),
                 d = as.integer(d),
                 K = as.integer(K),
                 mu = as.double(mu),
                 lambda_min = as.double(1),
                 obj_val = as.double(1))

    expect_equal(comp_result$lambda_min,result$lambda_min)
    expect_equal(comp_result$objective_value,result$obj_val)
    })


# test_that("FORCE Main Loop",{
#     set.seed(12345)
#     K <- 5
#     d <- 20
#     dat <- generate_glatent_dc(K,d,d,3,4)
#     sh <- t(dat$X)%*%dat$X / d
#     gh <- gforce.Gamma(dat$X)
#     diff <- diag(gh) - sh
#     initial_mixing <- 2/d
#     km_res <- gforce.kmeans(-D_Kmeans,K)
#     km_res <- km_res$clusters
#     km_sol <- gforce.clust2mat(km_res)

#     ren_start_res <- renegar_start(diff,K,initial_mixing,km_sol)
#     X <- ren_start_res$X1
#     E <- ren_start_res$E

#     E_EVEV <- eigen(E)
#     V <- E_EVEV$vectors
#     D <- diag(E_EVEV$values)
#     E_sqrt <- V%*%(D^(0.5))%*%t(V)
#     ESI <- solve(E_sqrt)

#     mu <- 0.5*0.01/log(d)
#     comp_result <- smoothed_objective(X,E,ESI,mu)

#     result <- .C("primal_dual_adar_R",
#                 D = as.double(diff),
#                 sh = as.double(sh),
#                 E = as.double(E),
#                 ESI = as.double(ESI),
#                 X= as.double(X),
#                 d = as.integer(d),
#                 K = as.integer(K),
#                 verbosity = as.integer(1),
#                 kmeans_iter = as.integer(10),
#                 max_iter = as.integer(200),
#                 finish_pgd = as.integer(1),
#                 number_restarts = as.integer(1),
#                 restarts = as.integer(50),
#                 alpha = as.double(0.001),
#                 Z_T = numeric(d^2),
#                 Z_best = numeric(d^2),
#                 km_best = as.integer(1:d),
#                 km_best_time = as.integer(1),
#                 km_iter_best = as.integer(1),
#                 km_iter_total = as.integer(1),
#                 dc = as.integer(1),
#                 dc_time = as.integer(1),
#                 dc_grad_iter = as.integer(1),
#                 grad_iter_best = as.integer(1),
#                 grad_iter_best_time = as.integer(1))

#     expect_equal(1,1)
#     })

#' @useDynLib GFORCE test_clust_to_opt_val
test_that("K-Means Objective Value",{

    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    opt_val_r <- 0
    result <- .C(test_clust_to_opt_val,
                 D = as.double(sh),
                 d = as.integer(d),
                 K = as.integer(K),
                 clusters = as.integer(dat$group_assignments),
                 opt_val = as.double(opt_val_r))

    opt_val <- sum(sh*gforce.clust2mat(dat$group_assignments))
    expect_equal(opt_val,result$opt_val)
    })
