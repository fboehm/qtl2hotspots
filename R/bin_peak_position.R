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

#' Binarize bin number
#'
#' @param bin_number a positive integer
#' @param total_number_of_bins a positive integer, at least as big as bin_number
#' @return a binary vector with exactly one nonzero entry
#' @export
#' @examples
#' binarize_bin_number(3, 10)
#' binarize_bin_number(10L, 10L)
#'
binarize_bin_number <- function(bin_number, total_number_of_bins){
  if (bin_number > total_number_of_bins) stop("bin_number must be less than or equal to bin_number")
  if (bin_number < 1) stop("bin_number must be at least 1")
  if (total_number_of_bins < 1) stop("total_number_of_bins must be at least 1")
  if (!is.numeric(bin_number)) stop("bin_number must be an integer")
  if (!is.numeric(total_number_of_bins)) stop("total_number_of_bins must be an integer")
  if (length(bin_number) != 1) stop("bin_number must have length 1")
  if (length(total_number_of_bins) != 1) stop("total_number_of_bins must have length 1")
  # end checks
  out <- rep(0, total_number_of_bins)
  out[bin_number] <- 1
  return(out)
}
