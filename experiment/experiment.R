rm(list=ls())
load("../results/factorization_exponential_families.RData")

i <- 1
table(sapply(res[[i]], length))
len_vec <- sapply(res[[i]], length)
res[[i]] <- res[[i]][which(len_vec > 1)]

###################

jj <- 1
x <- res[[i]][[jj]]
scalar_vec <- alpha_vec; family_val <- "curved_gaussian"
for(j in 1:length(x$fit)){
  nat_mat <- x$fit[[j]]$u_mat %*% t(x$fit[[j]]$v_mat)
  plot_prediction_against_observed(dat = x$dat, nat_mat_list = list(nat_mat),
                                         family = family_val, missing_idx_list = list(x$missing_idx),
                                         scalar = scalar_vec[j], plot = T)

}

# dat = x$dat
# nat_mat_list = list(nat_mat)
# width = 0.8
# family = family_val
# missing_idx_list = list(x$missing_idx)
# scalar = scalar_vec[j]
# plot = T
# max_points = 500000
#
# pred_mat_list <- lapply(nat_mat_list, function(nat_mat){
#   compute_mean(nat_mat, family = family, scalar = scalar)
# })
#
# tmp_list <- lapply(1:length(nat_mat_list), function(i){
#   cbind(dat[missing_idx_list[[i]]], pred_mat_list[[i]][missing_idx_list[[i]]])
# })
#
# angle_vec <- sapply(tmp_list, .compute_principal_angle)
# angle_val <- mean(angle_vec)
#
# tmp_mat <- do.call(rbind, tmp_list)
# if(nrow(tmp_mat) > max_points){
#   tmp_mat <- tmp_mat[sample(1:nrow(tmp_mat), max_points),]
# }
#
# res <- .within_prediction_region(tmp_mat, family = family, width = width, scalar = scalar, angle_val = angle_val)
# .plot_pca_diagnostic(tmp_mat, seq_vec = res$seq_vec, interval_mat = res$interval_mat,
#                      principal_line = res$principal_line, angle_val = angle_val)


###################

tmp_vec <- sapply(res[[i]], function(x){
  scalar_vec <- alpha_vec; family_val <- "curved_gaussian"

  quality_list <- lapply(1:length(x$fit), function(j){
    nat_mat <- x$fit[[j]]$u_mat %*% t(x$fit[[j]]$v_mat)
    eSVD::plot_prediction_against_observed(dat = x$dat, nat_mat_list = list(nat_mat),
                                           family = family_val, missing_idx_list = list(x$missing_idx),
                                           scalar = scalar_vec[j], plot = F)
  })

  quality_mat <- do.call(rbind, lapply(quality_list, unlist))
  if(any(quality_mat[,"bool"] == 1)){
    qualified_idx <- c(1:length(x$fit))[which(quality_mat[,"bool"] == 1)]
    qualified_idx[which.min(abs(quality_mat[qualified_idx, "angle_val"] - 45))]
  } else {
    which.min(abs(quality_mat[qualified_idx, "angle_val"] - 45))
  }
})

table(factor(as.numeric(tmp_vec), levels = c(1,2,3)))
