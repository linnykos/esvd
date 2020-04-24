set.seed(10)
load(paste0("../results/step0_screening", suffix, ".RData"))

max_val <- 5000
n <- nrow(dat); d <- ncol(dat)
cv_trials <- 2
svd_missing_list <- vector("list", length(cv_trials))

# generate list of missing indices
missing_idx_list <- lapply(1:cv_trials, function(j){
  set.seed(10*j)
  eSVD::construct_missing_values(n = n, p = d, num_val = 4)
})
training_idx_list <- lapply(1:cv_trials, function(j){
  c(1:prod(dim(dat)))[-missing_idx_list[[j]]]
})

log_dat <- log2(dat+1)

# fit many softImputes
starting_lambda <- min(sapply(1:cv_trials, function(i){
  log_dat_NA <- log_dat
  log_dat_NA[missing_idx_list[[i]]] <- NA

  softImpute::lambda0(log_dat_NA)
}))

paramMat_svd <- as.matrix(expand.grid(seq(1, starting_lambda, length.out = 50), c(5,10,20,50)))
colnames(paramMat_svd) <- c("lambda", "k")

svd_angle_res <- sapply(1:nrow(paramMat_svd), function(j){
  tmp_mat <- do.call(rbind, lapply(1:cv_trials, function(i){
    log_dat_NA <- log_dat
    log_dat_NA[missing_idx_list[[i]]] <- NA

    set.seed(10)
    softImpute_embedding <- softImpute::softImpute(log_dat_NA, rank.max = paramMat_svd[j, "k"], lambda = paramMat_svd[j, "lambda"])
    softImpute_pred <- softImpute_embedding$u %*% diag(softImpute_embedding$d) %*% t(softImpute_embedding$v)

    cbind(log_dat, softImpute_pred)
  }))

  training_val <- eSVD:::.compute_principal_angle(tmp_mat[missing_idx_list[[j]],])
  testing_val <- eSVD:::.compute_principal_angle(tmp_mat[missing_idx_list[[j]],])

  c(training = training_val, testing = testing_val)
})
svd_angle_res <- cbind(t(svd_angle_res), paramMat_svd)

# select the best model
idx <- which.min(abs(svd_angle_res[,2] - 45))
svd_missing <- lapply(1:cv_trials , function(i){
  log_dat_NA <- log_dat
  log_dat_NA[missing_idx_list[[i]]] <- NA

  set.seed(10)
  softImpute::softImpute(log_dat_NA, rank.max = paramMat_svd[idx, "k"], lambda = paramMat_svd[idx, "lambda"])
})


# fit the actual embedding as well
tmp <- svd(log_dat)
svd_embedding <- (n/d)^(1/4)*tmp$u[,1:paramMat_svd[idx, "k"]] %*% diag(sqrt(tmp$d[1:paramMat_svd[idx, "k"]]))

rm(list = c("init", "fit", "dat_NA", "j", "tmp", "starting_lambda", "idx",
            "log_dat", "log_dat_NA"))
print(paste0(Sys.time(), ": Finished naive SVD"))
save.image(paste0("../results/step1_naive_svd", suffix, ".RData"))

######################

# nat_mat_list <- lapply(1:length(svd_missing), function(i){
#   svd_missing[[i]]$u %*% diag(svd_missing[[i]]$d) %*% t(svd_missing[[i]]$v)
# })
#
# tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
#   cbind(dat[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
# }))
#
# pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
# ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))
# ang * 180/pi
