#' Convert a tibble with one list-column of vectors (of same length) to a tibble without list-column, where each entry of vector (within list-column) is a separate column.
#'
#' @param tibble a tibble containing one list-column
#' @param list_column name of list-column, enclosed with dplyr::quo()
#' @param another_column name of another column, enclosed with dplyr::quo()
#' @return a tibble with the list-column replaced with multiple columns. New columns have backticked numbers as names, like `1`
#' @examples
#' tibble(a = 1:3, b = c(list(1:2), list(3:4), list(5:6))) %>%
#' listcol_to_cols(list_column = dplyr::quo(b), another_column = dplyr::quo(a))
#' @export
listcol_to_cols <- function(tibble, list_column, another_column){
  tibble %>%
    tidyr::unnest() %>%
    dplyr::group_by(!! another_column) %>%
    dplyr::mutate(col=seq_along(!! another_column)) %>%
    tidyr::spread(key = col, value = !! list_column)
}
# see: https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
# see: https://stackoverflow.com/questions/50881440/split-a-list-column-into-multiple-columns-in-r
