#' Assign a peak position to a bin number
#'
#' @param peak_position one peak position (numeric)
#' @param breaks vector defining the boundaries of bins
#' @return an integer, the bin number
#' @export
#' @examples
#' bin_peak_position(5, 1:10)
#' bin_peak_position(1, 2)
#'
bin_peak_position <- function(peak_position, breaks){
  if (length(peak_position) > 1) stop("peak_position must be length 1")
  if (!is.numeric(peak_position)) stop("peak_position must be numeric")
  if (length(breaks) < 1) stop("breaks must have length at least 1")
  if (!is.numeric(breaks)) stop("breaks must be numeric")
  return(sum(peak_position > breaks) + 1)
}
