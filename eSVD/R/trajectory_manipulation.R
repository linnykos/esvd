.extract_pseudotimes <- function(slingshot_res){
  n <- length(slingshot_res$curves)

  mat_list <- lapply(1:n, function(curve_idx){
    # which ghost points are in the curve
    idx_all <- slingshot_res$idx[slingshot_res$curves[[curve_idx]]$idx]

    # all the pseudotimes
    lambda_all <- slingshot_res$curves[[curve_idx]]$lambda
    stopifnot(length(idx_all) == length(lambda_all))

    # which cells do these correspond to
    idx_cell <- sort(unique(idx_all))
    order_vec <- rep(NA, length(idx_cell))
    lambda_vec <- rep(NA, length(idx_cell))

    for(i in 1:length(idx_cell)){
      tmp_idx <- which(idx_all == idx_cell[i]) # idx of the ghost points
      order_vec[i] <- mean(which(slingshot_res$curves[[1]]$ord %in% tmp_idx))
      lambda_vec[i] <- mean(lambda_all[tmp_idx])
    }

    mat <- cbind(idx_cell,lambda_vec)
    colnames(mat) <- c("cell_idx", "pseudotime")
    mat
  })
}

.select_ideal_cells <- function(mat_list){
  stopifnot(length(mat_list) == 2)

  idx <- intersect(mat_list[[1]][,1], mat_list[[2]][,1])
  lambda_mat <- matrix(NA, nrow = length(idx), ncol = 2)
  for(i in 1:nrow(lambda_mat)){
    idx1 <- which(mat_list[[1]][,1] == idx[i])
    idx2 <- which(mat_list[[2]][,1] == idx[i])
    lambda_mat[i,] <- c(mat_list[[1]][idx1, 2], mat_list[[2]][idx2, 2])
  }

  max_val <- max(lambda_mat)
  acceptable_range <- max_val/10

  selected_idx <- which(abs(lambda_mat[,1] - lambda_mat[,2]) <= acceptable_range)
  selected_mat <- matrix(NA, nrow = length(selected_idx), ncol = 2)
  colnames(selected_mat) <- c("cell_idx", "pseudotime")
  selected_mat[,1] <- idx[selected_idx]
  selected_mat[,2] <- rowMeans(lambda_mat[selected_idx,])

  selected_mat
}

.append_trajectory_specific_cells <- function(pseudotime_mat, cell_mat, remaining_idx){
  stopifnot(remaining_idx %in% pseudotime_mat[,1])

  new_pseudotime_vec <- sapply(remaining_idx, function(i){
    pseudotime_mat[which(pseudotime_mat[,1] == i),2]
  })

  rbind(cell_mat, cbind(remaining_idx, new_pseudotime_vec))
}
