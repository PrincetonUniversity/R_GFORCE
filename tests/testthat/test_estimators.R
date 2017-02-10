context('Test PECOK Estimators')

# TEST FORM - test_that("NAME",{  })

#' @useDynLib GFORCE test_smoothed_gradient_S_base
test_that("Test V Measures",{
  K <- 5
  d <- 20
  dat <- generate_glatent_dc(K,d,d,3,4)
  gh1 <- gamma_hat(dat$X,par=TRUE)
  gh2 <- gamma_hat(dat$X,par=FALSE)

  expect_equal(gh1,gh2)
})