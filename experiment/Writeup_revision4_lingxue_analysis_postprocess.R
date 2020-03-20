rm(list=ls())
load("../results/lingxue_analysis.RData")

#fit_all_list <- fit_all_list[1:6]
lapply(fit_all_list, function(x){x$neg_bin_param})
lapply(fit_all_list, function(x){x$curved_gaussian_param})

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

###########################

for(i in 1:num_datasets){
  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup4_lingxue_embedding_", i, ".png"),
      height = 3000, width = 4000, res = 300,
      units = "px")

  par(mfrow = c(3,4))
  for(j in 1:3){
    main_vec <- c(paste0("Gaussian\n(constant variance, k = ", j+2, ")"),
                  paste0("Poisson,\n(k = ", j+2, ")"),
                  paste0("Negative binomial\n(k = ", j+2, ", est. size = ", fit_all_list[[i]][[j]][[3]]$neg_binom_param, ")"),
                  paste0("Curved Gaussian\n(k = ", j+2, ", est. alpha = ", fit_all_list[[i]][[j]][[4]]$curved_gaussian_param, ")"))

    for(k in 1:4){
      plot(fit_all_list[[i]][[j]][[k]]$u_mat[,1], fit_all_list[[i]][[j]][[k]]$u_mat[,2],
           asp = T, col = as.numeric(preprocessing_list[[i]]$label_vec), pch = 16,
           main = main_vec[k])
    }
  }

  graphics.off()
}

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
    scalar_vec <- c(sd(dat[fit_all_list[[i]][[j]]$missing_idx] - (fit_all_list[[i]][[j]][[1]]$u_mat %*% t(fit_all_list[[i]][[j]][[1]]$v_mat))[fit_all_list[[i]][[j]]$missing_idx]),
                    NA,
                    fit_all_list[[i]][[j]][[3]]$neg_binom_param,
                    fit_all_list[[i]][[j]][[4]]$curved_gaussian_param)
    main_vec <- c(paste0("Gaussian\n(constant variance, k = ", j+2, ")"),
                  paste0("Poisson,\n(k = ", j+2, ")"),
                  paste0("Negative binomial\n(k = ", j+2, ", est. size = ", fit_all_list[[i]][[j]][[3]]$neg_binom_param, ")"),
                  paste0("Curved Gaussian\n(k = ", j+2, ", est. alpha = ", fit_all_list[[i]][[j]][[4]]$curved_gaussian_param, ")"))


    for(k in 1:4){
      scalar <- scalar_vec[k]
      nat_mat <- fit_all_list[[i]][[j]][[k]]$u_mat %*% t(fit_all_list[[i]][[j]][[k]]$v_mat)

      plot_prediction_against_observed(dat = dat, nat_mat = nat_mat,
                                       family = distribution_vec[k], scalar = scalar,
                                       missing_idx = fit_all_list[[i]][[j]]$missing_idx,
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
    scalar_vec <- c(sd(dat[fit_all_list[[i]][[j]]$missing_idx] - (fit_all_list[[i]][[j]][[1]]$u_mat %*% t(fit_all_list[[i]][[j]][[1]]$v_mat))[fit_all_list[[i]][[j]]$missing_idx]),
                    NA,
                    fit_all_list[[i]][[j]][[3]]$neg_binom_param,
                    fit_all_list[[i]][[j]][[4]]$curved_gaussian_param)
    main_vec <- c(paste0("Gaussian\n(constant variance, k = ", j+2, ")"),
                  paste0("Poisson,\n(k = ", j+2, ")"),
                  paste0("Negative binomial\n(k = ", j+2, ", est. size = ", fit_all_list[[i]][[j]][[3]]$neg_binom_param, ")"),
                  paste0("Curved Gaussian\n(k = ", j+2, ", est. alpha = ", fit_all_list[[i]][[j]][[4]]$curved_gaussian_param, ")"))


    for(k in 1:4){
      scalar <- scalar_vec[k]
      nat_mat <- fit_all_list[[i]][[j]][[k]]$u_mat %*% t(fit_all_list[[i]][[j]][[k]]$v_mat)

      plot_prediction_against_observed(dat = dat, nat_mat = nat_mat,
                                       family = distribution_vec[k], scalar = scalar,
                                       missing_idx = c(1:prod(dim(dat)))[-fit_all_list[[i]][[j]]$missing_idx],
                                       main = main_vec[k])
    }
  }
  graphics.off()
}






