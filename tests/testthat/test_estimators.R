context('Test PECOK Estimators')

# TEST FORM - test_that("NAME",{  })

test_that("Test V Measures",{
  K <- 5
  d <- 20
  dat <- generate_glatent_dc(K,d,d,3,4)
  gh1 <- gforce.Gamma(dat$X,par=TRUE)
  gh2 <- gforce.Gamma(dat$X,par=FALSE)

  expect_equal(gh1,gh2)
})