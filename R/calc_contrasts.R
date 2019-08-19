#' Calculate all contrasts for 8 founder alleles
#'
#' @param grouping a numeric vector, length 8, assigning founder alleles to each of 3 groups: 1, 2, or 3.
#' @param effects founder allele effects vector of length 8
#'
#' @export
calc_contrasts <- function(grouping, effects){
  means <- numeric(length = 2) # one mean per founder allele
  for (i in 1:2){
    means[i] <- mean(effects_vector[which(grouping == i)])
    names(means)[i] <- LETTERS[which(grouping == i)]
  }
  diff <- abs(means[1] - means[2])
  return(tibble::tibble(contrast = contr, value = diff))
}
