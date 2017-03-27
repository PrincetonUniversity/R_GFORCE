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
  expect_equal(0,max(max(abs(fs2-t(fs2))))) #symmetry
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
  expect_equal(0,max(max(abs(fs1-t(fs1))))) #symmetry
  expect_true(min(eigen(fs1)$values) > 0) # Positive Definite
  expect_true(min(eigen(fs2)$values) > 0) # Positive Definite

  })

test_that("Renegar Method Initial Solution -- R Implementation",{
  d <- 50
  K <- 7
  s <- 0.25
  dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
  D <- - t(dat$X)%*%dat$X / d
  opt_estimate <- kmeans(D,K)$cluster
  opt_estimate <- gforce.clust2mat(opt_estimate)
  rs <- gforce.FORCE.init(D,K,s,opt_estimate,R_only=TRUE)
  
  # test strict feasibility of E
  E <- rs$E
  expect_equal(0,max(abs(colSums(E) - 1))) #E 1 = 1
  expect_equal(K,sum(diag(E))) #trace == K
  expect_equal(0,max(max(abs(E-t(E))))) #symmetry
  expect_true(min(eigen(E)$values) > 0) # Positive Definite

  # test feasibility of X
  X0 <- rs$X0
  expect_equal(0,max(abs(colSums(X0) - 1))) #E 1 = 1
  expect_equal(K,sum(diag(X0))) #trace == K
  expect_equal(0,max(max(abs(X0-t(X0))))) #symmetry
  expect_true(min(eigen(X0)$values) > 0) # Positive Definite

  # test objective values
  X0_val <- sum(D*X0)
  E_val <- sum(D*E)
  expect_equal(rs$E_val,E_val)
  expect_equal(rs$X0_val,X0_val)
  expect_true(X0_val < E_val)
  })

test_that("Renegar Method Initial Solution -- C Implementation",{
  d <- 50
  K <- 7
  s <- 0.25
  dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
  D <- - t(dat$X)%*%dat$X / d
  opt_estimate <- kmeans(D,K)$cluster
  opt_estimate <- gforce.clust2mat(opt_estimate)
  rs <- gforce.FORCE.init(D,K,s,opt_estimate,R_only=FALSE)
  
  # test strict feasibility of E
  E <- rs$E
  expect_equal(0,max(abs(colSums(E) - 1))) #E 1 = 1
  expect_equal(K,sum(diag(E))) #trace == K
  expect_equal(0,max(max(abs(E-t(E))))) #symmetry
  expect_true(min(eigen(E)$values) > 0) # Positive Definite

  # test feasibility of X
  X0 <- rs$X0
  expect_equal(0,max(abs(colSums(X0) - 1))) #E 1 = 1
  expect_equal(K,sum(diag(X0))) #trace == K
  expect_equal(0,max(max(abs(X0-t(X0))))) #symmetry
  expect_true(min(eigen(X0)$values) > 0) # Positive Definite

  # test objective values
  X0_val <- sum(D*X0)
  E_val <- sum(D*E)
  expect_equal(rs$E_val,E_val)
  expect_equal(rs$X0_val,X0_val)
  expect_true(X0_val < E_val)
  })

test_that("Renegar Method Initial Solution -- Cluster Representation",{
  d <- 50
  K <- 7
  s <- 0.25
  dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
  D <- - t(dat$X)%*%dat$X / d
  rs <- gforce.FORCE.init(D,K,s,dat$group_assignments,cluster_representation=TRUE,R_only=FALSE)
  
  # test strict feasibility of E
  E <- rs$E
  expect_equal(0,max(abs(colSums(E) - 1))) #E 1 = 1
  expect_equal(K,sum(diag(E))) #trace == K
  expect_equal(0,max(max(abs(E-t(E))))) #symmetry
  expect_true(min(eigen(E)$values) > 0) # Positive Definite

  # test feasibility of X
  X0 <- rs$X0
  expect_equal(0,max(abs(colSums(X0) - 1))) #E 1 = 1
  expect_equal(K,sum(diag(X0))) #trace == K
  expect_equal(0,max(max(abs(X0-t(X0))))) #symmetry
  expect_true(min(eigen(X0)$values) > 0) # Positive Definite

  # test objective values
  X0_val <- sum(D*X0)
  E_val <- sum(D*E)
  expect_equal(rs$E_val,E_val)
  expect_equal(rs$X0_val,X0_val)
  expect_true(X0_val < E_val)
  })

# test_that("Renegar Method Initial Solution -- Cluster Representation",{
#   d <- 50
#   K <- 7
#   s <- 0.25
#   dat <- gforce.generator(K,d,d,3,graph='DeltaC',cov_gap_mult=4)
#   D <- - t(dat$X)%*%dat$X / d
#   rs <- gforce.FORCE.init(D,K,s,dat$group_assignments,cluster_representation=TRUE,R_only=TRUE)
  
#   # test strict feasibility of E
#   E <- rs$E
#   expect_equal(0,max(abs(colSums(E) - 1))) #E 1 = 1
#   expect_equal(K,sum(diag(E))) #trace == K
#   expect_equal(0,max(max(abs(E-t(E))))) #symmetry
#   expect_true(min(eigen(E)$values) > 0) # Positive Definite

#   # test feasibility of X
#   X0 <- rs$X0
#   expect_equal(0,max(abs(colSums(X0) - 1))) #E 1 = 1
#   expect_equal(K,sum(diag(X0))) #trace == K
#   expect_equal(0,max(max(abs(X0-t(X0))))) #symmetry
#   expect_true(min(eigen(X0)$values) > 0) # Positive Definite

#   # test objective values
#   X0_val <- sum(D*X0)
#   E_val <- sum(D*E)
#   expect_equal(rs$E_val,E_val)
#   expect_equal(rs$X0_val,X0_val)
#   expect_true(X0_val < E_val)
#   })


#' @useDynLib GFORCE test_random_shuffle
test_that("Shuffle Indices",{
  n <-1000
  shuffle_sum <- sum(0:(n-1))
  result <- .C(test_random_shuffle,
                n=as.integer(n),
                shuffled=as.integer(numeric(n)))
  expect_equal(sum(result$shuffled),shuffle_sum)
  })