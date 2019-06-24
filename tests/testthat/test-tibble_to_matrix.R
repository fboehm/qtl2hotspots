library(hotspots)
library(testthat)

foo <- tibble::tibble(v1 = letters[c(1, 1, 2, 2)], v2 = letters[c(3, 4, 3, 4)], lrt = runif(4))

test_that("tibble converts to matrix with correct dimensions", {
  expect_equal(nrow(tibble_to_matrix(foo)), 2)
  expect_equal(ncol(tibble_to_matrix(foo)), 2)
})
