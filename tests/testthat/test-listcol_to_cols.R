library(qtl2hotspots)
library(testthat)

myfoo <- tibble::tibble(num = 1:4, numlist = list(1:10, 11:20, 21:30, 31:40))
myout <- listcol_to_cols(tibble = myfoo,
                         list_column = dplyr::quo(numlist),
                         another_column = dplyr::quo(num)
                         )


test_that("correct dimensions & rownames in simple cases", {
  expect_equal(dim(myout), c(4, 11))
  expect_equal(colnames(myout)[1], "num")
  expect_equal(colnames(myout)[2:11], as.character(1:10))
})

test_that("non-quosures give errors", {
  expect_error(listcol_to_cols(myfoo, num, numlist))
})
