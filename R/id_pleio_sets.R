#' Cluster pleiotropy test statistics
#'
#' @param tib a tibble with three columns; the first two are gene ids and the last is pleiotropy test statistics
#' @param binarize logical of length one indicating whether to binarize the matrix of pleiotropy test statistics based on pleio_threshold value
#' @param pleio_threshold a positive number to distinguish pleiotropy from separate QTL

cluster_pleio_stats <- function(tib, binarize = FALSE, pleio_threshold = NULL){
  # convert tibble to matrix
  mat <- tibble_to_matrix(tib)
  # binarize matrix, optionally
  if (binarize) mat <- mat < pleio_threshold
  # sort matrix
  hclust(dist(mat))$order -> ord
  bmat[ord, ord]
}


#' Convert tibble to matrix
#'
#' @param tib

tibble_to_matrix <- function(tib){
  rn <- unique(tib[[1]])
  cn <- unique(tib[[2]])
  out <- matrix(NA, nrow = length(rn), ncol = length(cn))
  rownames(out) <- rn
  colnames(out) <- cn
  for (i in 1:length(rn)){
    for (j in 1:length(cn)){
      out[i, j] <- tib[[3]][(tib[[1]] == rownames(out)[i] & tib[[2]] == colnames(out)[j])]
    }
  }
  return(out)
}

#' Sort a binary matrix into something close to a block diagonal
#'
#' @param bmat a binary matrix

sort_binary_matrix <- function(bmat){
  bmat[ , order(colSums(bmat), decreasing = TRUE)]
}
