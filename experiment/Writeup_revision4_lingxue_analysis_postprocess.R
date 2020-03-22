rm(list=ls())
load("../results/lingxue_analysis.RData")

num_datasets <- 6
neg_binom_selected <- vector("list", num_datasets)
curved_gaussian_selected <- vector("list", num_datasets)

# select best neg binom or curved gaussian
for(i in 1:num_datasets){ # loop over datasets
  tmp_vec <- rep(NA, 3)

  for(k in 1:3){
    vec <- sapply(1:7, function(j){
      dat_impute <- preprocessing_list[[i]]$dat_impute

      nat_mat_list <- lapply(1:2, function(x){
        u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[k]]$neg_binom_missing[[j]]$u_mat
        v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[k]]$neg_binom_missing[[j]]$v_mat
        u_mat  %*% t(v_mat)
      })
      missing_idx_list <- lapply(1:2, function(x){
        fit_all_list[[(x-1)*num_datasets+i]][[k]]$missing_idx
      })

      plot_prediction_against_observed(dat_impute, nat_mat_list, family = "neg_binom", missing_idx_list,
                                       scalar = neg_binom_vec[j], plot = F)
    })

    tmp_vec[k] <- which.min(abs(vec - 45))
  }

  neg_binom_selected[[i]] <- tmp_vec
}

for(i in 1:num_datasets){ # loop over datasets
  tmp_vec <- rep(NA, 3)

  for(k in 1:3){
    vec <- sapply(1:7, function(j){
      dat_impute <- preprocessing_list[[i]]$dat_impute

      nat_mat_list <- lapply(1:2, function(x){
        u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[k]]$curved_gaussian_missing[[j]]$u_mat
        v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[k]]$curved_gaussian_missing[[j]]$v_mat
        u_mat  %*% t(v_mat)
      })
      missing_idx_list <- lapply(1:2, function(x){
        fit_all_list[[(x-1)*num_datasets+i]][[k]]$missing_idx
      })

      plot_prediction_against_observed(dat_impute, nat_mat_list, family = "curved_gaussian", missing_idx_list,
                                       scalar = curved_gaussian_vec[j], plot = F)
    })

    tmp_vec[k] <- which.min(abs(vec - 45))
  }

  curved_gaussian_selected[[i]] <- tmp_vec
}

# ###########################
#
# for(i in 1:num_datasets){
#   png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_lingxue_embedding_", i, ".png"),
#       height = 3000, width = 4000, res = 300,
#       units = "px")
#
#   par(mfrow = c(3,4))
#   for(j in 1:3){
#     main_vec <- c(paste0("Gaussian\n(constant variance, k = ", j+2, ")"),
#                   paste0("Poisson,\n(k = ", j+2, ")"),
#                   paste0("Negative binomial\n(k = ", j+2, ", est. size = ", fit_all_list[[i]][[j]][[3]]$neg_binom_param, ")"),
#                   paste0("Curved Gaussian\n(k = ", j+2, ", est. alpha = ", fit_all_list[[i]][[j]][[4]]$curved_gaussian_param, ")"))
#
#     for(k in 1:4){
#       plot(fit_all_list[[i]][[j]][[k]]$u_mat[,1], fit_all_list[[i]][[j]][[k]]$u_mat[,2],
#            asp = T, col = as.numeric(preprocessing_list[[i]]$label_vec), pch = 16,
#            main = main_vec[k])
#     }
#   }
#
#   graphics.off()
# }

##########
distribution_vec <- c("gaussian", "poisson", "neg_binom", "curved_gaussian")

# testing
for(i in 1:num_datasets){
  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_lingxue_test_", i, ".png"),
      height = 3000, width = 4000, res = 300,
      units = "px")
  dat <- preprocessing_list[[i]]$dat_impute

  par(mfrow = c(3,4))
  for(j in 1:3){
    tmp_mat <- rbind(cbind(dat[fit_all_list[[i]][[j]]$missing_idx],
                           (fit_all_list[[i]][[j]][[1]]$u_mat %*% t(fit_all_list[[i]][[j]][[1]]$v_mat))[fit_all_list[[i]][[j]]$missing_idx]),
                     cbind(dat[fit_all_list[[i+6]][[j]]$missing_idx],
                           (fit_all_list[[i+6]][[j]][[1]]$u_mat %*% t(fit_all_list[[i+6]][[j]][[1]]$v_mat))[fit_all_list[[i+6]][[j]]$missing_idx]))

    scalar_vec <- c(sd(tmp_mat[,1] - tmp_mat[,2]),
                    NA,
                    neg_binom_vec[neg_binom_selected[[i]][j]],
                    curved_gaussian_vec[curved_gaussian_selected[[i]][j]])
    main_vec <- c(paste0("Gaussian\n(constant variance, k = ", j+2, ")"),
                  paste0("Poisson,\n(k = ", j+2, ")"),
                  paste0("Negative binomial\n(k = ", j+2, ", est. size = ", scalar_vec[3], ")"),
                  paste0("Curved Gaussian\n(k = ", j+2, ", est. alpha = ", scalar_vec[4], ")"))


    for(k in 1:4){
      scalar <- scalar_vec[k]

      if(k %in% c(1,2)){
        nat_mat_list <- lapply(1:2, function(x){
          u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]][[k]]$u_mat
          v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]][[k]]$v_mat
          u_mat  %*% t(v_mat)
        })
      } else if(k == 3){
        nat_mat_list <- lapply(1:2, function(x){
          u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$neg_binom_missing[[neg_binom_selected[[i]][j]]]$u_mat
          v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$neg_binom_missing[[neg_binom_selected[[i]][j]]]$v_mat
          u_mat  %*% t(v_mat)
        })
      } else {
        nat_mat_list <- lapply(1:2, function(x){
          u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$curved_gaussian_missing[[curved_gaussian_selected[[i]][j]]]$u_mat
          v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$curved_gaussian_missing[[curved_gaussian_selected[[i]][j]]]$v_mat
          u_mat  %*% t(v_mat)
        })
      }

      missing_idx_list <- lapply(1:2, function(x){
        fit_all_list[[(x-1)*num_datasets+i]][[j]]$missing_idx
      })

      plot_prediction_against_observed(dat = dat, nat_mat_list = nat_mat_list,
                                       family = distribution_vec[k], scalar = scalar,
                                       missing_idx_list = missing_idx_list,
                                       main = main_vec[k])
    }
  }
  graphics.off()
}


# training
for(i in 1:num_datasets){
  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_lingxue_train_", i, ".png"),
      height = 3000, width = 4000, res = 300,
      units = "px")
  par(mfrow = c(3,4))
  dat <- preprocessing_list[[i]]$dat_impute

  for(j in 1:3){
    tmp_mat <- rbind(cbind(dat[fit_all_list[[i]][[j]]$missing_idx],
                           (fit_all_list[[i]][[j]][[1]]$u_mat %*% t(fit_all_list[[i]][[j]][[1]]$v_mat))[fit_all_list[[i]][[j]]$missing_idx]),
                     cbind(dat[fit_all_list[[i+6]][[j]]$missing_idx],
                           (fit_all_list[[i+6]][[j]][[1]]$u_mat %*% t(fit_all_list[[i+6]][[j]][[1]]$v_mat))[fit_all_list[[i+6]][[j]]$missing_idx]))

    scalar_vec <- c(sd(tmp_mat[,1] - tmp_mat[,2]),
                    NA,
                    neg_binom_vec[neg_binom_selected[[i]][j]],
                    curved_gaussian_vec[curved_gaussian_selected[[i]][j]])
    main_vec <- c(paste0("Gaussian\n(constant variance, k = ", j+2, ")"),
                  paste0("Poisson,\n(k = ", j+2, ")"),
                  paste0("Negative binomial\n(k = ", j+2, ", est. size = ", scalar_vec[3], ")"),
                  paste0("Curved Gaussian\n(k = ", j+2, ", est. alpha = ", scalar_vec[4], ")"))


    for(k in 1:4){
      scalar <- scalar_vec[k]

      if(k %in% c(1,2)){
        nat_mat_list <- lapply(1:2, function(x){
          u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]][[k]]$u_mat
          v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]][[k]]$v_mat
          u_mat  %*% t(v_mat)
        })
      } else if(k == 3){
        nat_mat_list <- lapply(1:2, function(x){
          u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$neg_binom_missing[[neg_binom_selected[[i]][j]]]$u_mat
          v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$neg_binom_missing[[neg_binom_selected[[i]][j]]]$v_mat
          u_mat  %*% t(v_mat)
        })
      } else {
        nat_mat_list <- lapply(1:2, function(x){
          u_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$curved_gaussian_missing[[curved_gaussian_selected[[i]][j]]]$u_mat
          v_mat <- fit_all_list[[(x-1)*num_datasets+i]][[j]]$curved_gaussian_missing[[curved_gaussian_selected[[i]][j]]]$v_mat
          u_mat  %*% t(v_mat)
        })
      }

      missing_idx_list <- lapply(1:2, function(x){
        c(1:prod(dim(dat)))[-fit_all_list[[(x-1)*num_datasets+i]][[j]]$missing_idx]
      })

      plot_prediction_against_observed(dat = dat, nat_mat_list = nat_mat_list,
                                       family = distribution_vec[k], scalar = scalar,
                                       missing_idx_list = missing_idx_list,
                                       main = main_vec[k])
    }
  }
  graphics.off()
}






