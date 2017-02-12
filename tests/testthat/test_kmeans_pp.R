
context('Test K-means++ Algorithm')

# TEST FORM - test_that("NAME",{  })

#' @useDynLib GFORCE test_smoothed_gradient_S_base
test_that("GS_t Base",{

    K <- 5
    d <- 20
    dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
    sh <- t(dat$X)%*%dat$X / d
    initial_mixing <- 2/d
    km_res <- kmeanspp(-sh,K)

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

#' @useDynLib GFORCE kmeans_pp_R
test_that("K-means++",{

    K <- 5
    d <- 20
    num_random_exps <- 100
    count_not_perfect <- 0
    total_purity <- 0

    for(i in 1:num_random_exps){
        dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
        sh <- t(dat$X)%*%dat$X / 20
        result <- .C(kmeans_pp_R,
                     D=as.double(sh),
                     K= as.integer(K),
                     n = as.integer(d),
                     m = as.integer(d),
                     group_assignments = as.integer(rep(0,d)),
                     centers = numeric(K*d))
        purity_experiment <- purity_measure(result$group_assignments,dat$group_assignments)
        if(purity_experiment != 1){
            count_not_perfect <- count_not_perfect + 1
        }
        total_purity <- total_purity + purity_experiment
    }
    # cat(sprintf('K-means++ Test: %.2f average purity\r\n',total_purity/num_random_exps))
    # cat(sprintf('K-means++ Test: %d / %d perfect recovery\r\n',num_random_exps-count_not_perfect,num_random_exps))
    expect_false(is.nan(total_purity))
    })