library(hotspots)
library(testthat)


test_that("bin assignment is correct in simple cases", {
  expect_equal(bin_peak_position(1, 2), 1)
  expect_equal(bin_peak_position(2, 1), 2)
  expect_equal(bin_peak_position(5, 1:10), 5)
})


test_that("errors occur when input is bad", {
  expect_error(bin_peak_position(1:2, 2))
  expect_error(bin_peak_position("a", 2))
  expect_error(bin_peak_position(1, NULL))
  expect_error(bin_peak_position(1, "a"))
})
