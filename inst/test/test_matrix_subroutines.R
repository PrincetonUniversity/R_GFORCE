context('Test BLAS/LAPACK-Alike Matrix Routines')

# TEST FORM - test_that("NAME",{  })
test_that("daps",{
    d <- 20
    A <- runif(d)
    c <- runif(1)

    Apc_comp <- A + c

    result <- .C("test_daps",
                 A = as.double(A),
                 c = as.double(c),
                 d = as.integer(d))

    expect_equal(result$A,Apc_comp)
    })

test_that("dsmtd",{
    d <- 20
    A <- matrix(runif(d^2),ncol=d)
    B <- runif(d)

    AB <- A %*% diag(B)
    BA <- diag(B) %*% A

    resultL <- .C("test_dsmtd",
                 A = as.double(A),
                 B = as.double(B),
                 d = as.integer(d),
                 side = as.character('L'))
    C_BA <- matrix(resultL$A,ncol=d)

    resultR <- .C("test_dsmtd",
                 A = as.double(A),
                 B = as.double(B),
                 d = as.integer(d),
                 side = as.character('R'))
    C_AB <- matrix(resultR$A,ncol=d)
    
    expect_equal(C_BA,BA)
    expect_equal(C_AB,AB)
    })

test_that("dvexp",{
    d <- 100
    A <- runif(d)
    A_exp_comp <- exp(A)

    result <- .C("test_dvexp",
                 A = as.double(A),
                 d = as.integer(d))
    A_exp <- result$A

    expect_equal(A_exp_comp,A_exp)
    })


test_that("dsumv",{
    d <- 100
    A <- runif(d)
    sum_r <- 0

    result <- .C("test_dsumv",
                 A = as.double(A),
                 d = as.integer(d),
                 sum_r = as.double(sum_r))
    
    expect_equal(sum(A),result$sum_r)
    })

test_that("dtrace",{
    d <- 20
    A <- matrix(runif(d^2),ncol=d)
    trace_r <- 0

    result <- .C("test_dtrace",
                 A = as.double(A),
                 d = as.integer(d),
                 trace_r = as.double(trace_r))
    
    expect_equal(sum(diag(A)),result$trace_r)
    })


test_that("dcsum",{
    d <- 20
    A <- matrix(runif(d^2),ncol=d)

    result <- .C("test_dcsum",
                 A = as.double(A),
                 d = as.integer(d),
                 A_csums = numeric(d))
    
    expect_equal(colSums(A),result$A_csums)
    })

test_that("dxpyez",{
    d <- 20
    A <- matrix(runif(d^2),ncol=d)
    B <- matrix(runif(d^2),ncol=d)

    result <- .C("test_dxpyez",
                 d = as.integer(d^2),
                 A = as.double(A),
                 B = as.double(B),
                 C = numeric(d^2))
    C <- matrix(result$C,ncol=d)
    expect_equal(C,A + B)
    })

