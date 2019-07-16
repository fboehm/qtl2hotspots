#' If the inputted number exceeds `threshold`, return `threshold`; otherwise, return the inputted value.
#'
#' @param input numeric vector to be modified
#' @param threshold threshold value
#' @export
#' @return a numeric vector with length equaling that of `input`
#' @examples
#' trunc2(input = 1:10, threshold = 5)
trunc2 <- function(input, threshold){
  input[input > threshold] <- threshold
  return(input)
}
