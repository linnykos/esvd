#' Construct indices for missing values
#'
#' @param n number of rows of intended matrix
#' @param p number of columns of intended matrix
#' @param num_val number of values to make missing in each row/column
#'
#' @return vector of indices between 1 and \code{n*p}
#' @export
construct_missing_values <- function(n, p, num_val = ceiling(min(c(4, n/10, p/10)))){
  missing_mat <- rbind(do.call(rbind, (lapply(1:n, function(x){
    cbind(x, sample(1:p, num_val))
  }))), do.call(rbind, (lapply(1:p, function(x){
    cbind(sample(1:n, num_val), x)
  }))))

  unique(apply(missing_mat, 1, function(x){
    (x[2]-1)*n + x[1]
  }))
}
