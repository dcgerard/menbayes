test_that("unknown parental genotypes works", {
  expect_no_error({
  polymapr_approx_g(x = c(5, 5, 0, 0, 0), g1 = NULL, g2 = NULL)
  polymapr_approx_g(x = c(0, 5, 0, 10, 0), g1 = NULL, g2 = NULL)
  })
})
