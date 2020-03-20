rm(list=ls())
load("../results/lingxue_analysis.RData")

# interactive 3D plot
# https://rpubs.com/aagarwal29/179912

num_datasets <- 7

# select best neg binom or curved gaussian
for(i in 1:num_datasets){ # loop over datasets
  for(k in 1:3){
    vec <- sapply(1:7, function(j){
      dat_impute <- preprocessing_list[[i]]$dat_impute
      u_mat <- fit_all_list[[i]][[k]]$neg_binom_missing[[j]]$u_mat
      v_mat <- fit_all_list[[i]][[k]]$neg_binom_missing[[j]]$v_mat
      nat_mat <- u_mat  %*% t(v_mat)
      missing_idx <- fit_all_list[[i]][[k]]$missing_idx

      plot_prediction_against_observed(dat_impute, nat_mat, family = "neg_binom", missing_idx = missing_idx,
                                       scalar = neg_binom_vec[j], plot = F)
    })

    idx <- which.min(abs(vec - 45))
    tmp <- fit_all_list[[i]][[k]]$neg_binom_missing[[idx]]
    tmp$neg_binom_param <- neg_binom_vec[idx]
    fit_all_list[[i]][[k]][[3]] <- tmp
  }
}

for(i in 1:num_datasets){ # loop over datasets
  for(k in 1:3){
    vec <- sapply(1:7, function(j){
      dat_impute <- preprocessing_list[[i]]$dat_impute
      u_mat <- fit_all_list[[i]][[k]]$curved_gaussian_missing[[j]]$u_mat
      v_mat <- fit_all_list[[i]][[k]]$curved_gaussian_missing[[j]]$v_mat
      nat_mat <- u_mat  %*% t(v_mat)
      missing_idx <- fit_all_list[[i]][[k]]$missing_idx

      plot_prediction_against_observed(dat_impute, nat_mat, family = "curved_gaussian", missing_idx = missing_idx,
                                       scalar = curved_gaussian_vec[j], plot = F)
    })

    idx <- which.min(abs(vec - 45))
    tmp <- fit_all_list[[i]][[k]]$curved_gaussian_missing[[idx]]
    tmp$curved_gaussian_param <- curved_gaussian_vec[idx]
    fit_all_list[[i]][[k]][[4]] <- tmp
  }
}

#############

i <- 4 # which dataset?
j <- 1 # how many dimensions?
k <- 3 # which distribution?

u_mat <- fit_all_list[[i]][[j]][[k]]$u_mat

rgl::plot3d(u_mat[,1], u_mat[,2], u_mat[,3], asp = T, col = as.numeric(preprocessing_list[[i]]$label_vec),
            main = "Dataset 4")
