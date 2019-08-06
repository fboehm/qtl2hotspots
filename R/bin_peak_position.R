#' Assign a peak position to a bin number
#'
#' @param peak_position one peak position (numeric)
#' @param breaks vector defining the boundaries of bins
#' @return an integer, the bin number
#' @export
#'
bin_peak_position <- function(peak_position, breaks){
  return(sum(peak_position > breaks) + 1)
}
