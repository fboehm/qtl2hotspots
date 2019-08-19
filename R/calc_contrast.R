#' Calculate all contrasts for 8 founder alleles
#'
#' @param grouping a numeric vector, length 8, assigning founder alleles to each of 3 groups: 1, 2, or 3.
#' @param effects founder allele effects vector of length 8
#'
#' @export
calc_contrast <- function(grouping, effects){
  # return NA if contrast can't be calculated, say, for a vector of all 1's.
  if (length(unique(grouping[grouping < 3])) <= 1) {
    out <- NA
    } else {
      means <- numeric(length = 2) # one mean per founder allele
      for (i in 1:2){
        means[i] <- mean(effects_vector[which(grouping == i)])
        names(means)[i] <- paste(LETTERS[which(grouping == i)], collapse = "")
        }
      diff <- abs(means[1] - means[2])
      contr <- paste(names(means), collapse = "_")
      out <- tibble::tibble(contrast = contr, value = diff)
    }
  return(out)
}
