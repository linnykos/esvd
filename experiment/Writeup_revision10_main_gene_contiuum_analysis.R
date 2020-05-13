rm(list=ls())
# load("../results/tmp.RData")
#
# cluster_labels <- as.numeric(cell_type_vec)
# order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
# cluster_group_list <- lapply(order_vec, function(x){
#   grep(paste0("^", x), levels(cell_type_vec))
# })
#
# upscale_factor <- 1
# p <- 3
#
# set.seed(10)
# esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                                cluster_group_list = cluster_group_list,
#                                verbose = T, upscale_factor = upscale_factor,
#                                reduction_percentage = 0.2,
#                                squared = T)
# esvd_curves$lineages
#
# ######################################################
#
# # focus on one trajectory for now
# curve_idx <- 2
# idx_all <- esvd_curves$idx[esvd_curves$curves[[curve_idx]]$idx] # which ghost points are in the first lineage
# lambda_all <- esvd_curves$curves[[curve_idx]]$lambda
# stopifnot(length(idx_all) == length(lambda_all))
# idx_cell <- sort(unique(idx_all)) # which cells do these correspond to
# order_vec <- rep(NA, length(idx_cell))
# lambda_vec <- rep(NA, length(idx_cell))
#
# for(i in 1:length(idx_cell)){
#   tmp_idx <- which(idx_all == idx_cell[i]) # idx of the ghost points
#   order_vec[i] <- mean(which(esvd_curves$curves[[1]]$ord %in% tmp_idx))
#   lambda_vec[i] <- mean(lambda_all[tmp_idx])
# }
#
# nat_mat <- esvd_embedding$u_mat %*% t(esvd_embedding$v_mat)
# pred_mat <- eSVD::compute_mean(nat_mat = nat_mat, family = "curved_gaussian")
#
# circular_list <- vector("list", length = ncol(pred_mat))
#
# for(j in 1:ncol(pred_mat)){
#   print(j)
#   if(j %% floor(ncol(pred_mat)/10) == 0) {
#     save.image("../results/tmp_continuum.RData")
#   }
#
#   vec <- pred_mat[idx_cell[order(order_vec)], j]
#   circular_list[[j]] <- eSVD::find_highly_expressed_region(vec, resolution = 50)
# }
#
# save.image("../results/tmp_continuum.RData")

load("../results/tmp_continuum.RData")
circular_list2 <- vector("list", length = ncol(pred_mat))
for(j in 1:ncol(pred_mat)){
  print(j)
  if(j %% floor(ncol(pred_mat)/10) == 0) {
    save.image("../results/tmp_continuum2.RData")
  }

  vec <- dat_impute[idx_cell[order(order_vec)], j]
  circular_list[[j]] <- eSVD::find_highly_expressed_region(vec, resolution = 50)
}
save.image("../results/tmp_continuum2.RData")

