library(hotspots)
library(testthat)

foo <- tibble::tibble(v1 = letters[c(1, 1, 2, 2)], v2 = letters[c(3, 4, 3, 4)], lrt = runif(4))

test_that("tibble converts to rectangular matrix with correct dimensions", {
  expect_equal(nrow(tibble_to_matrix(foo, symmetric = FALSE)), 2)
  expect_equal(ncol(tibble_to_matrix(foo, symmetric = FALSE)), 2)
})

bar <- foo
bar$lrt[3] <- NA


test_that("tibble converts to symmetric matrix with correct dimensions", {
  expect_equal(nrow(tibble_to_matrix(bar, symmetric = TRUE)), 4)
  expect_equal(ncol(tibble_to_matrix(foo, symmetric = TRUE)), 4)
})

test_that("tibble converts to symmetric matrix with correct dimensions", {
  expect_equal(sum(is.na(tibble_to_matrix(bar, symmetric = TRUE))), 10)
  expect_equal(sum(is.na(tibble_to_matrix(foo, symmetric = TRUE))), 8)
  expect_true(isSymmetric(tibble_to_matrix(foo, symmetric = TRUE)))
})
