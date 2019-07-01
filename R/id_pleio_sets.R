#' Cluster pleiotropy test statistics
#'
#' @param tib a tibble with three columns; the first two are gene ids and the last is pleiotropy test statistics
#' @param binarize logical of length one indicating whether to binarize the matrix of pleiotropy test statistics based on pleio_threshold value
#' @param pleio_threshold a positive number to distinguish pleiotropy from separate QTL
#' @export

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
#' @param tib a tibble with pleiotropy test results
#' @param symmetric logical to indicate if the resulting matrix should be symmetric
#' @export

tibble_to_matrix <- function(tib, symmetric = FALSE){
  if (!symmetric){
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
  }
  if (symmetric) {
    rn <- unique(union(tib[[1]], tib[[2]]))
    out <- matrix(NA, nrow = length(rn), ncol = length(rn))
    rownames(out) <- rn
    colnames(out) <- rn
    for (i in 2:length(rn)){
      for (j in 1:(i - 1)){
        ind1 <- tib[[1]] == rownames(out)[i] & tib[[2]] == colnames(out)[j]
        ind2 <- tib[[2]] == rownames(out)[i] & tib[[1]] == colnames(out)[j]
        if (sum(ind1) | sum(ind2)) out[i, j] <- tib[[3]][ind1 | ind2]
        out[j, i] <- out[i, j]
      }
    }
  }
  return(out)
}

#' Determine adjacency matrix from a n local by m nonlocal traits binary matrix
#'
#' @param mat a binary matrix with rownames and colnames
#' @export

make_adjacency_matrix <- function(mat){
  mat %*% t(mat)
}
