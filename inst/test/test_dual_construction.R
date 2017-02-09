

context('Test Dual Construction')

# TEST FORM - test_that("NAME",{  })
test_that("Dual Construction",{
    set.seed(12345) #need to set seed bc doesnt always exist
    dat <- generate_glatent_dc(5,20,20,3,4)
    gh <- gamma_hat(dat$X)
    sh <- t(dat$X)%*%dat$X / 20
    diff <- sh - diag(gh)
    diff <- -1*diff

    K <- 5
    d <- 20
    eps1 <- 0.01
    eps2 <- 10^(-7)
    Y_T_min <- 0.01
    #Y_a <- rep(0,d)
    Y_T <- 0
    feasible <- 0

    result <- .C("kmeans_dual_solution_primal_min_R",
                 ga_hat=as.integer(dat$group_assignments),
                 D=as.double(diff),
                 K= as.integer(K),
                 dimension=as.integer(d),
                 eps1 =as.double(eps1),
                 eps2 =as.double(eps2),
                 Y_T_min =as.double(Y_T_min),
                 Y_a = numeric(d),
                 Y_T =as.double(Y_T),
                 feasible = as.integer(feasible))

    dc <- dual_solution(dat$group_assignments,-diff)
    # print(result$Y_a)
    # print(result$Y_T)

    expect_equal(sum(result$Y_a)*2 + K*result$Y_T,dc$dual_value)
    })