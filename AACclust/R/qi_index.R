
#' Quality Index (QI) for Partition Assessment
#'
#' Computes the QI Index that evaluates partition quality through pairwise
#' comparison of observations.
#'
#' @param data A data.frame or matrix of categorical variables.
#' @param grupos A list of integer vectors representing the partition.
#'
#' @return A numeric value in [0, 1] representing the partition quality.
#'
#' @examples
#' data(felinos)
#' partition <- aac_cluster(felinos, alpha = 5)
#' qi_index(felinos, partition)
#'
#' @export
qi_index <- function(data, grupos) {
  if (length(grupos) == 0) return(NA)
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  n <- ncol(data_mat)
  
  C <- calcular_matriz_C_cpp(data_mat)
  
  S <- matrix(0, m, m)
  for (k in seq_along(grupos)) {
    idx <- grupos[[k]]
    S[idx, idx] <- 1
  }
  
  termino1 <- sum((3 * C - 2 * n) * S)
  termino2 <- sum(2 * (n - C))
  numerador <- termino1 + termino2
  denominador <- (2 - 1/m) * m^2 * n
  
  return(numerador / denominador)
}

