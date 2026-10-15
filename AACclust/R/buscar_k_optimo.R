
#' Search for Optimal Number of Clusters
#'
#' Evaluates the QI Index for all k in [1, m] using parallel computing.
#'
#' @param data A data.frame or matrix of categorical variables.
#' @param n_cores Number of cores to use. Defaults to detectCores() - 1.
#'
#' @return A list with:
#'   \item{k_optimo}{Optimal number of clusters}
#'   \item{QI_maximo}{Maximum QI value}
#'   \item{QI_por_k}{Vector of QI values for each k}
#'   \item{particiones}{List of partitions for each k}
#'
#' @examples
#' data(felinos)
#' result <- buscar_k_optimo(felinos)
#' print(result$k_optimo)
#'
#' @export
#' @importFrom parallel mclapply detectCores
Find_k_optimo <- function(data, n_cores = NULL) {
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  n <- ncol(data_mat)
  
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  C_global <- calcular_matriz_C_cpp(data_mat)
  
  resultados <- parallel::mclapply(1:m, function(k) {
    if (k == 1) {
      grupos <- list(1:m)
    } else if (k == m) {
      grupos <- as.list(1:m)
    } else {
      grupos <- alfa19_cpp(data_mat, as.integer(k))
    }
    qi <- calcular_QI_particion_cpp(C_global, grupos, m, n)
    list(k = k, QI = qi, grupos = grupos)
  }, mc.cores = n_cores)
  
  QI_por_k <- sapply(resultados, function(x) x$QI)
  particiones <- lapply(resultados, function(x) x$grupos)
  
  k_optimo <- which.max(QI_por_k)
  
  list(
    k_optimo = k_optimo,
    QI_maximo = QI_por_k[k_optimo],
    QI_por_k = QI_por_k,
    particiones = particiones
  )
}

