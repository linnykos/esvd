.find_neighbor <- function(mat, i, num_neighbors){
  row_vec <- c(1:nrow(mat))[-i]
  idx_me <- which(mat[i,] == 0)
  vec <- sapply(row_vec, function(x){
    idx_other <- which(mat[x,] == 0)
    length(intersect(idx_me, idx_other))/length(unique(c(idx_me, idx_other)))
  })

  row_vec[order(vec, decreasing = T)[1:num_neighbors]]
}

.find_true_zeros <- function(dropout_mat, num_neighbors = NA){
  if(is.na(num_neighbors)) num_neighbors <- ceiling(nrow(dropout_mat)/10)

  neighbor_list <- lapply(1:nrow(dropout_mat), function(x){
    .find_neighbor(dropout_mat, x, num_neighbors)
  })

  zero_mat <- dropout_mat
  for(i in 1:nrow(dropout_mat)){
    idx <- which(dropout_mat[i,] == 0)
    vec <- rep(0, length(idx))
    zz <- colSums(dropout_mat[neighbor_list[[i]], idx])
    vec[which(zz != 0)] <- NA
    zero_mat[i, idx] <- vec
  }

  zero_mat
}
