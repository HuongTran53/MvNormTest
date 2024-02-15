test_that("Skewness estimator works", {
  set.seed(1);
  x <- rnorm(100)
  expect_equal(round(skewness(x), 3), -0.0720)
})
