rm(list=ls())
load("../results/step5_trajectory_original.RData")

num_order_vec_esvd <- c(5, rep(3,2), c(6,1,1,1,6,4), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec2_esvd <- color_func(0.3)[num_order_vec_esvd]
col_name_esvd <- c("orange", rep("bluish green", 2), c("gray", "yellow", "yellow", "yellow", "gray", "blue"), rep("skyblue", 2), rep("orange", 2))
order_vec_esvd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_esvd <- data.frame(name = levels(cell_type_vec),
                            idx = sort(unique(cluster_labels)),
                            order = order_vec_esvd,
                            col_name = col_name_esvd,
                            col_code = col_vec_esvd)
col_info_esvd$factor_idx <- as.numeric(as.factor(col_info_esvd$col_name))
col_info_esvd[,c(5,6)] <- col_info_esvd[,c(6,5)]
colnames(col_info_esvd)[c(5,6)] <- colnames(col_info_esvd)[c(6,5)]
col_info_esvd

###########

num_clust <- c(3,3,3,3, 3,3,3,3, 3,3,3,3, 3)

# # cluster each cluster label
# i <- 13
# k_vec <- c(2:20)
# obj_vec <- rep(NA, length(k_vec))
# idx <- which(cluster_labels == i)
# for(k in 1:length(k_vec)){
#   set.seed(10)
#   res <- stats::kmeans(esvd_embedding$u_mat[idx,], centers = k)
#   obj_vec[k] <- res$tot.withinss
# }
# plot(k_vec, obj_vec, pch = 16)

# construct cell centers across all cell types
cluster_center_list <- vector("list", length(num_clust))
for(i in 1:max(cluster_labels)){
  idx <- which(cluster_labels == i)
  set.seed(10)
  res <- stats::kmeans(esvd_embedding$u_mat[idx,], centers = num_clust[i])
  cluster_center_list[[i]] <- res$centers
}

combn_mat <- combn(3,2)
par(mfrow = c(1,3))
for(i in 1:ncol(combn_mat)){
  dim1 <- combn_mat[1,i]; dim2 <- combn_mat[2,i]
  plot(esvd_embedding$u_mat[,dim1], esvd_embedding$u_mat[,dim2], asp = T, pch = 16,
       col = col_vec2_esvd[cluster_labels])

  for(k in 1:length(cluster_center_list)){
    cluster_centers <- cluster_center_list[[k]]
    points(cluster_centers[,dim1], cluster_centers[,dim2], pch = 16, cex = 3, col = "black")
    points(cluster_centers[,dim1], cluster_centers[,dim2], pch = 16, cex = 2, col = col_info_esvd$col_code[k])
  }
}

# so now, compute the inner product between each gene's latent vector and each cluster center
average_expression_list <- vector("list", length(cluster_center_list))
for(i in 1:length(cluster_center_list)){
  average_expression_list[[i]] <- esvd_embedding$v_mat %*% t(cluster_center_list[[i]])
  average_expression_list[[i]] <- 1/(average_expression_list[[i]])
}

# standardize all the values
for(j in 1:ncol(dat_impute)){
  vec_all <- unlist(lapply(average_expression_list, function(mat){mat[j,]}))
  sd_val <- stats::sd(vec_all)
  mean_val <- mean(vec_all)
  for(k in 1:length(average_expression_list)){
    average_expression_list[[k]][j,] <- (average_expression_list[[k]][j,]-mean_val)/sd_val
  }
}

# find highly expressed values
obj_vec_list <- vector("list", length(cluster_center_list))
for(k in 1:length(cluster_center_list)){
  max_vec <- apply(average_expression_list[[k]], 1, max)
  other_max <- apply(sapply(average_expression_list[-k], function(mat){apply(mat, 1, max)}), 1, max)

  obj_vec_list[[k]] <- max_vec-other_max
}
sapply(obj_vec_list, max)




