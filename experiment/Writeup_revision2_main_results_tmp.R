rm(list=ls())
load("../results/step3_scalar_heuristic_tmp.RData")

# missing_idx_list <- lapply(1:cv_trials, function(j){
#   set.seed(10*j)
#   missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
#     cbind(x, sample(1:d, 4))
#   }))), do.call(rbind, (lapply(1:d, function(x){
#     cbind(sample(1:n, 4), d)
#   }))))
#
#   dat_impute_NA <- dat_impute
#   for(tmp in 1:nrow(missing_idx)){
#     dat_impute_NA[missing_idx[tmp,1], missing_idx[tmp,2]] <- NA
#   }
#
#   which(is.na(dat_impute_NA))
# })

sapply(res_list, length)
res_list <- res_list[-9]

for(i in 1:length(res_list)){
  mat_list <- lapply(1:length(res_list[[i]]), function(j){
    # pred_mat <- res_list[[i]][[j]]$u_mat %*% t(res_list[[i]][[j]]$v_mat)
    pred_mat <- 1/(res_list[[i]][[j]]$u_mat %*% t(res_list[[i]][[j]]$v_mat))

    cbind(dat_impute[missing_idx_list[[j]]], pred_mat[missing_idx_list[[j]]])
    # cbind(log(dat_impute[missing_idx_list[[j]]]+1), log(pred_mat[missing_idx_list[[j]]]+1))
  })

  mat <- do.call(rbind, mat_list)

  pca_res <- stats::prcomp(mat, center = F, scale = F)
  diag_vec <- c(1,1)

  plot(mat[,2], mat[,1], asp = T, pch = 16, col = rgb(0,0,0,0.2))
  lines(c(0,1e6), c(0,1e6), col = "red", lwd = 2)
  lines(c(0, 1e6), c(0, 1e6*pca_res$rotation[1,1]/pca_res$rotation[2,1]), col = "blue", lwd = 2, lty = 2)

  x_val <- seq(1, 1e5, length.out = 100)
  y_val_top <- sapply(x_val, function(x){stats::qnorm(0.75, mean = x, sd = x/scalar_vec[i])})
  y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.25, mean = x, sd = x/scalar_vec[i])})

  polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
          border = NA, density = 30, angle = -45)
  lines(x_val, y_val_top, col = "red", lwd = 2, lty = 2)
  lines(x_val, y_val_bottom, col = "red", lwd = 2, lty = 2)
}

plot(res_list[[3]][[1]]$u_mat[,1], res_list[[3]][[1]]$u_mat[,2], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
plot(res_list[[3]][[1]]$u_mat[,1], res_list[[3]][[1]]$u_mat[,3], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
plot(res_list[[3]][[1]]$u_mat[,2], res_list[[3]][[1]]$u_mat[,3], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
