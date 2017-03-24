context('Test Initial Solution Construction')

test_that("Test Full Rank Feasible R vs Native",{
    K <- 3
    d <- 20
    fs1 <- gforce.full_rank_feasible(d,K,R_only=FALSE)
    fs2 <- gforce.full_rank_feasible(d,K,R_only=TRUE)
    expect_equal(fs1,fs2)

    K <- 4
    fs1 <- gforce.full_rank_feasible(d,K,R_only=FALSE)
    fs2 <- gforce.full_rank_feasible(d,K,R_only=TRUE)
    expect_equal(fs1,fs2)

    })

test_that("Test Feasibility",{
    K <- 3
    d <- 20
    fs1 <- gforce.full_rank_feasible(d,K,R_only=FALSE)
    fs2 <- gforce.full_rank_feasible(d,K,R_only=TRUE)
    expect_equal(0,max(abs(colSums(fs1) - 1))) #E 1 = 1
    expect_equal(0,max(abs(colSums(fs2) - 1))) #E 1 = 1
    expect_equal(K,sum(diag(fs1))) #trace == K
    expect_equal(K,sum(diag(fs2))) #trace == K
    expect_equal(0,max(max(abs(fs1-t(fs1))))) #symmetry
    expect_true(min(eigen(fs1)$values) > 0) # Positive Definite
    expect_true(min(eigen(fs2)$values) > 0) # Positive Definite

    K <- 4
    fs1 <- gforce.full_rank_feasible(d,K,R_only=FALSE)
    fs2 <- gforce.full_rank_feasible(d,K,R_only=TRUE)
    expect_equal(0,max(abs(colSums(fs1) - 1))) #E 1 = 1
    expect_equal(0,max(abs(colSums(fs2) - 1))) #E 1 = 1
    expect_equal(K,sum(diag(fs1))) #trace == K
    expect_equal(K,sum(diag(fs2))) #trace == K
    expect_equal(0,max(max(abs(fs1-t(fs1))))) #symmetry
    expect_true(min(eigen(fs1)$values) > 0) # Positive Definite
    expect_true(min(eigen(fs2)$values) > 0) # Positive Definite

    })

