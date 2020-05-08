rm(list=ls())
load("../results/step2_baron_scalar_tuning.RData")
paramMat_esvd <- paramMat; colnames(paramMat_esvd)[3] <- "scalar"

for(dat_i in 1:length(preprocessing_list)){
  print(dat_i)

  svd_missing_list <- svd_missing_list_list[[dat_i]]
  dat <- preprocessing_list[[dat_i]]$dat_impute
  dim(dat)
  missing_idx_list <- missing_idx_list_list[[dat_i]]

  nat_mat_list <- lapply(1:length(svd_missing_list), function(i){
    svd_missing_list[[i]]$u %*% diag(svd_missing_list[[i]]$d) %*% t(svd_missing_list[[i]]$v)
  })

  tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
    cbind(log2(dat + 1)[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
  }))
  sd_val <- sd(tmp_mat[,1] - tmp_mat[,2])

  training_idx_list <- lapply(1:length(missing_idx_list), function(i){
    c(1:prod(dim(dat)))[-missing_idx_list[[i]]]
  })

  png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision10_baron_training_testing_", dat_i, ".png"),
      height = 1500, width = 2500, res = 300,
      units = "px")
  par(mfrow = c(1,2))
  plot_prediction_against_observed(log2(dat+1), nat_mat_list = nat_mat_list,
                                   missing_idx_list = training_idx_list,
                                   family = "gaussian", scalar = sd_val,
                                   main = "SVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                   max_points = 1e6)


  plot_prediction_against_observed(log2(dat+1), nat_mat_list = nat_mat_list,
                                   missing_idx_list = missing_idx_list,
                                   family = "gaussian", scalar = sd_val,
                                   main = "SVD embedding:\nMatrix-completion diagnostic\n(Testing set)")
  graphics.off()
}

############################


# for each latent dimension and distr, pick the appropriate scalar
for(dat_i in 1:5){
  print(paste0("Working on dataset ", dat_i))
  dat <- preprocessing_list[[dat_i]]$dat_impute
  dat <- dat * 1000/max(dat)
  missing_idx_list <- missing_idx_list_list[[dat_i]]
  esvd_missing_list <- esvd_missing_list_list[[dat_i]]

  nat_mat_list_list <- lapply(1:nrow(paramMat_esvd), function(i){
    lapply(1:cv_trials, function(j){
      u_mat <- esvd_missing_list[[i]][[j]]$u_mat
      v_mat <- esvd_missing_list[[i]][[j]]$v_mat
      u_mat %*% t(v_mat)
    })
  })

  training_idx_list <- lapply(1:length(missing_idx_list), function(i){
    c(1:prod(dim(dat)))[-missing_idx_list[[i]]]
  })

  selection_mat <- as.matrix(expand.grid(k_vec, c(2,3)))
  tuning_idx <- rep(NA, nrow(selection_mat))
  for(j in 1:nrow(selection_mat)){
    distr_num <- selection_mat[j,2]
    k <- selection_mat[j,1]
    param_idx <- intersect(which(paramMat_esvd[,"fitting_distr"] == distr_num), which(paramMat_esvd[,"k"] == k))

    fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[distr_num]
    esvd_angle_res <- eSVD::tuning_select_scalar(dat = dat, nat_mat_list_list = nat_mat_list_list[param_idx],
                                                 family = fitting_distr,  missing_idx_list = missing_idx_list,
                                                 scalar_vec = paramMat_esvd[param_idx,"scalar"])
    tuning_idx[j] <- esvd_angle_res$idx
  }

  # now plot
  png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision10_baron_testing_esvd_", dat_i, ".png"),
      height = 2000, width = 3000, res = 300,
      units = "px")
  par(mfrow = c(2,3))
  for(j in 1:nrow(selection_mat)){
    print(j)

    distr_num <- selection_mat[j,2]
    k <- selection_mat[j,1]
    param_idx <- intersect(which(paramMat_esvd[,"fitting_distr"] == distr_num), which(paramMat_esvd[,"k"] == k))
    fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[distr_num]
    if(fitting_distr == "neg_binom"){
      scalar <- neg_binom_vec[tuning_idx[j]]
    } else {
      scalar <- curved_gaussian_vec[tuning_idx[j]]
    }
    main_str <- paste0(ifelse(fitting_distr == "neg_binom", "Negative binomial", "Curved Gaussian"),
                       "\n(k = ", k, ", scalar = ", scalar, ")")

    plot_prediction_against_observed(dat, nat_mat_list = nat_mat_list_list[param_idx][[tuning_idx[j]]],
                                     missing_idx_list = missing_idx_list,
                                     family = fitting_distr, scalar = scalar,
                                     main = main_str,
                                     max_points = 1e6)
  }
  graphics.off()
}

