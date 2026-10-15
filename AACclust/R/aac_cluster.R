
#' Affinity Contraction Agglomerative Clustering (AAC)
#'
#' Clustering algorithm for categorical data based on exact attribute matching
#' and mutual affinity edge contractions.
#'
#' @param data A data.frame or matrix of categorical variables.
#' @param alpha Integer. Target number of clusters.
#'
#' @return A list of integer vectors, each containing the indices of
#'   observations belonging to that cluster.
#'
#' @examples
#' data(felinos)
#' partition <- aac_cluster(felinos, alpha = 5)
#' print(partition)
#'
#' @export
#' @useDynLib AACclust, .registration = TRUE
#' @importFrom Rcpp sourceCpp
aac_cluster <- function(data, alpha) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data.frame or matrix")
  }
  if (!is.numeric(alpha) || alpha < 1) {
    stop("'alpha' must be a positive integer")
  }
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "double"
  alfa19_cpp(data_mat, as.integer(alpha))
}

