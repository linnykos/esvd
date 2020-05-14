construct_pseudotime_trajectory_matrix <- function(slingshot_res, cluster_labels){
  stopifnot(length(slingshot_res$curves) == 2) #only works currently for 2 branching lineages
  len <- length(slingshot_res$curves)

  # first extract each trajectories set of pseudotimes
  pseudotime_df_list <- .extract_pseudotimes(slingshot_res)

  # next, focus on cells shared between both trajectories
  shared_df <- .compile_common_cells(pseudotime_df_list)

  # deal with other cells
  specific_df <- .compile_unique_cells(pseudotime_df_list)

  # merge the data frames
  tmp <- cbind(do.call(rbind, specific_df), NA)
  colnames(tmp) <- colnames(shared_df)
  all_df <- rbind(shared_df, tmp)

  # all cluster labels
  all_df$cluster_labels <- cluster_labels[all_df$cell_idx]
  all_df <- all_df[order(all_df$cell_idx),]

  all_df
}

####################################

.extract_pseudotimes <- function(slingshot_res){
  n <- length(slingshot_res$curves)

  df_list <- lapply(1:n, function(curve_idx){
    # which ghost points are in the curve
    idx_all <- slingshot_res$idx[slingshot_res$curves[[curve_idx]]$idx]

    # all the pseudotimes
    lambda_all <- slingshot_res$curves[[curve_idx]]$lambda
    dist_all <- slingshot_res$curves[[curve_idx]]$dist_ind
    stopifnot(length(idx_all) == length(lambda_all))

    # which cells do these correspond to
    idx_cell <- sort(unique(idx_all))
    lambda_vec <- rep(NA, length(idx_cell))
    dist_vec <- rep(NA, length(idx_cell))

    for(i in 1:length(idx_cell)){
      tmp_idx <- which(idx_all == idx_cell[i]) # idx of the ghost points
      lambda_vec[i] <- mean(lambda_all[tmp_idx])
      dist_vec[i] <- mean(dist_all[tmp_idx])
    }

    data.frame(cell_idx = idx_cell, pseudotime = lambda_vec, dist_to_curve = dist_vec)
  })
}

.compile_common_cells <- function(df_list){
  stopifnot(length(df_list) == 2) # hard-coded

  idx <- intersect(df_list[[1]]$cell_idx, df_list[[2]]$cell_idx)
  tol_val <- max(c(df_list[[1]]$pseudotime, df_list[[2]]$pseudotime))/10
  lambda_vec <- rep(NA, length(idx))
  bool_vec <- rep(NA, length(idx))

  for(i in 1:length(idx)){
    idx1 <- which(df_list[[1]]$cell_idx == idx[i])
    idx2 <- which(df_list[[2]]$cell_idx == idx[i])

    lambda1 <- df_list[[1]]$pseudotime[idx1]
    lambda2 <- df_list[[2]]$pseudotime[idx2]

    if(df_list[[1]]$dist_to_curve[idx1] < df_list[[2]]$dist_to_curve[idx2]){
      lambda_vec[i] <- lambda1
    } else {
      lambda_vec[i] <- lambda2
    }

    bool_vec[i] <- (abs(lambda1 - lambda2) <= tol_val)
  }

  data.frame(cell_idx = idx, pseudotime = lambda_vec, consensus = bool_vec)
}

.compile_unique_cells <- function(pseudotime_df_list){
  stopifnot(length(pseudotime_df_list) == 2) # hard-coded
  len <- length(pseudotime_df_list)

  lapply(1:len, function(i){
    other_branch <- ifelse(i == 1, 2, 1)
    specific_idx <- which(!pseudotime_df_list[[i]]$cell_idx %in%
                                 pseudotime_df_list[[other_branch]]$cell_idx)
    pseudotime_df_list[[i]][specific_idx, c(1:2)]
  })
}
