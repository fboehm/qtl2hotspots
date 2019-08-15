library(qtl2hotspots)
library(testthat)


test_that("binary vector is correct in simple cases", {
  expect_equal(sum(binarize_bin_number(1, 2)), 1)
  expect_length(binarize_bin_number(1, 2), 2)
  expect_equal(class(binarize_bin_number(1, 2)), "numeric")
})

