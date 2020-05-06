set.seed(10)
load(paste0("../results/step0_baron_screening", suffix, ".RData"))

k_vec <- c(3,5,10)
cv_trials <- 3
missing_idx_list_list <- lapply(1:length(preprocessing_list), function(i){
  n <- nrow(preprocessing_list[i]$dat_impute)
  p <-  nrow(preprocessing_list[i]$dat_impute)

  missing_idx_list <- lapply(1:cv_trials, function(j){
    set.seed(10*j)
    eSVD::construct_missing_values(n = n, p = p, num_val = 2)
  })

  missing_idx_list
})

###############################

# determine the best tuning parameter for each dataset
svd_tuning_list <- vector("list", length(preprocessing_list))

for(i in 1:length(preprocessing_list)){
  print(paste0("Tuning dataset ", i))
  dat <- preprocessing_list[[i]]$dat_impute
  missing_idx_list <- missing_idx_list_list[[i]]

  training_idx_list <- lapply(1:cv_trials, function(j){
    c(1:prod(dim(dat)))[-missing_idx_list[[j]]]
  })

  log_dat <- log2(dat+1)

  # fit many softImputes
  starting_lambda <- min(sapply(1:cv_trials, function(j){
    log_dat_NA <- log_dat
    log_dat_NA[missing_idx_list[[j]]] <- NA

    softImpute::lambda0(log_dat_NA)
  }))

  paramMat_svd <- as.matrix(expand.grid(seq(1, starting_lambda, length.out = 25), k_vec))
  colnames(paramMat_svd) <- c("lambda", "k")

  svd_angle_res <- sapply(1:nrow(paramMat_svd), function(j){
    tmp_res <- sapply(1:cv_trials, function(ii){
      log_dat_NA <- log_dat
      log_dat_NA[missing_idx_list[[ii]]] <- NA

      set.seed(10)
      softImpute_embedding <- softImpute::softImpute(log_dat_NA, rank.max = paramMat_svd[j, "k"],
                                                     lambda = paramMat_svd[j, "lambda"])
      softImpute_pred <- softImpute_embedding$u %*% diag(softImpute_embedding$d) %*% t(softImpute_embedding$v)

      tmp_mat <- cbind(as.numeric(log_dat), as.numeric(softImpute_pred))

      training_val <- eSVD:::.compute_principal_angle(tmp_mat[training_idx_list[[ii]],])
      testing_val <- eSVD:::.compute_principal_angle(tmp_mat[missing_idx_list[[ii]],])

      c(training = training_val, testing = testing_val)
    })

    rowMeans(tmp_res)
  })

  svd_angle_res <- cbind(t(svd_angle_res), paramMat_svd)

  # select the best model
  idx <- which.min(abs(svd_angle_res[,2] - 45))

  lambda_val <- paramMat_svd[idx, "lambda"]
  k_val <- paramMat_svd[idx, "k"]

  svd_tuning_list[[i]] <- list(lambda = lambda_val, k = k_val, svd_angle_res = svd_angle_res)
  save.image(paste0("../results/step1_baron_gaussian_fitting", suffix, ".RData"))
}

print(paste0(Sys.time(), ": Finished selecting tuning parameters"))

#############################

# now fit the actual missing value SVDs
svd_missing_list_list <- vector("list", length(preprocessing_list))

for(i in 1:length(preprocessing_list)){
  print(paste0("Fitting missing value for dataset ", i))
  dat <- preprocessing_list[[i]]$dat_impute
  missing_idx_list <- missing_idx_list_list[[i]]

  log_dat <- log2(dat+1)

  tmp <- lapply(1:cv_trials , function(j){
    log_dat_NA <- log_dat
    log_dat_NA[missing_idx_list[[j]]] <- NA

    set.seed(10)
    softImpute::softImpute(log_dat_NA, rank.max = svd_tuning_list[[i]]$k, lambda = svd_tuning_list[[i]]$lambda)
  })

  svd_missing_list_list[[i]] <- tmp
}

print(paste0(Sys.time(), ": Finished fitting missing value SVD"))
save.image(paste0("../results/step1_baron_gaussian_fitting", suffix, ".RData"))

# now fit the actual SVDs
svd_embedding_list <- lapply(1:length(preprocessing_list), function(i){
  dat <- preprocessing_list[[i]]$dat_impute
  log_dat <- log2(dat+1)

  svd(log_dat)
})

print(paste0(Sys.time(), ": Finished naive SVD"))
save.image(paste0("../results/step1_baron_gaussian_fitting", suffix, ".RData"))
